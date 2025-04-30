from typing import List, Dict, Tuple, Set
import copy
import time
import math
import networkx as nx
from pysat.solvers import Glucose4
import dfg
import graph
import main  # Thêm import main để truy cập trực tiếp vào namespace của main

# Các import cũ giữ nguyên
from main import create_array, create_control_signal_array
from main import add_initial_conditions, add_final_conditions, add_existence_constraints
from main import add_communication_constraints, add_calculation_constraints
from main import add_block_constraints, add_control_signal_capacity, add_capacity_constraints
from main import add_clause, interpret_solution

# Khai báo biến toàn cục
DFG = None
CGRA = None

"""
Module này triển khai phương pháp phân cấp cho ánh xạ DFG lên CGRA.
"""

class HierarchicalMapper:
    def __init__(self, dfg_obj, cgra_obj, max_cycles=None, nregs=2, nopes=1, num_partitions=2):
        self.DFG = dfg_obj
        self.CGRA = cgra_obj
        self.nregs = nregs
        self.nopes = nopes
        self.max_cycles = max_cycles if max_cycles else dfg_obj.ndata * 2
        self.num_partitions = num_partitions
        self.partitions = []
        self.interfaces = []
        self.node_to_partition = {}
        self.local_mappings = []
        self.global_mapping = None
        self.solve_time = 0
        
        # Đảm bảo DFG có phương thức get_ndata
        if not hasattr(self.DFG, 'get_ndata'):
            self.DFG.get_ndata = lambda: self.DFG.ndata
            
        # Đảm bảo DFG có phương thức get_ninputs
        if not hasattr(self.DFG, 'get_ninputs'):
            self.DFG.get_ninputs = lambda: self.DFG.ninputs
            
        # Bổ sung thêm thuộc tính nexs nếu thiếu
        if not hasattr(self.DFG, 'nexs'):
            self.DFG.nexs = 0

    def partition_dfg(self) -> List[List[int]]:
        """
        Phân vùng DFG thành các phân vùng nhỏ hơn có kích thước tương đương
        và mối liên kết giữa các phân vùng ít nhất.
        """
        print(f"Partitioning DFG into {self.num_partitions} partitions...")
        start_time = time.time()
        
        # Tạo NetworkX graph từ DFG
        G = nx.DiGraph()
        for i in range(self.DFG.get_ndata()):
            G.add_node(i)
        
        # Thêm cạnh dựa trên operands
        for i in range(self.DFG.get_ndata()):
            for operand_set in self.DFG.operands[i]:
                for operand in operand_set:
                    if operand >= 0:  # Bỏ qua operand âm
                        G.add_edge(operand, i)
        
        # Tính toán mức của mỗi node
        levels = {}
        sources = [i for i in range(self.DFG.get_ninputs())]
        
        for source in sources:
            levels[source] = 0
        
        # BFS để tính mức
        remaining_nodes = list(range(self.DFG.get_ndata()))
        for source in sources:
            if source in remaining_nodes:
                remaining_nodes.remove(source)
        
        while remaining_nodes:
            for node in list(remaining_nodes):
                predecessors = list(G.predecessors(node))
                if all(pred not in remaining_nodes for pred in predecessors):
                    level = max([levels.get(pred, 0) for pred in predecessors], default=0) + 1
                    levels[node] = level
                    remaining_nodes.remove(node)
        
        # Nhóm node theo mức
        level_groups = {}
        for node, level in levels.items():
            if level not in level_groups:
                level_groups[level] = []
            level_groups[level].append(node)
        
        # Phân phối node vào các phân vùng
        partitions = [[] for _ in range(self.num_partitions)]
        sorted_levels = sorted(level_groups.keys())
        
        # Chia theo level để đảm bảo tính phụ thuộc
        nodes_per_partition = self.DFG.get_ndata() // self.num_partitions
        remainder = self.DFG.get_ndata() % self.num_partitions
        
        # Tính cân đối các phân vùng
        target_sizes = [nodes_per_partition + (1 if i < remainder else 0) for i in range(self.num_partitions)]
        
        # Phân phối các node theo level vào các phân vùng
        current_partition = 0
        for level in sorted_levels:
            node_group = level_groups[level]
            # Sắp xếp node theo độ kết nối (nodes với nhiều kết nối đặt cùng phân vùng)
            node_group.sort(key=lambda n: len(list(G.successors(n))) + len(list(G.predecessors(n))), reverse=True)
            
            for node in node_group:
                # Tìm partition phù hợp nhất cho node
                best_partition = self._find_best_partition(node, partitions, G, target_sizes)
                partitions[best_partition].append(node)
                self.node_to_partition[node] = best_partition
                
                # Cập nhật target size
                target_sizes[best_partition] -= 1
        
        # Kiểm tra và cân bằng lại nếu cần
        self._balance_partitions(partitions, G, target_sizes)
        
        self.partitions = partitions
        partition_time = time.time() - start_time
        print(f"Partitioning completed in {partition_time:.3f} seconds")
        
        # In thông tin các phân vùng
        for i, partition in enumerate(self.partitions):
            print(f"Partition {i}: {len(partition)} nodes")
        
        return self.partitions

    def _find_best_partition(self, node, partitions, graph, target_sizes):
        """
        Tìm phân vùng tốt nhất cho node dựa trên:
        1. Số lượng kết nối với các node trong phân vùng
        2. Kích thước hiện tại của phân vùng
        """
        scores = []
        
        # Tính điểm cho mỗi phân vùng
        for i, partition in enumerate(partitions):
            if target_sizes[i] <= 0:  # Phân vùng đã đầy
                scores.append(-1000)
                continue
                
            # Tính số lượng kết nối với partition
            connections = 0
            for p_node in partition:
                if graph.has_edge(node, p_node) or graph.has_edge(p_node, node):
                    connections += 1
            
            # Điểm = (số kết nối * 2 - số lượng node hiện tại)
            score = connections * 2 - len(partition)
            scores.append(score)
        
        # Nếu tất cả các partition có điểm âm (đã đầy), chọn partition ít node nhất
        if all(score < 0 for score in scores):
            return min(range(len(partitions)), key=lambda i: len(partitions[i]))
            
        return scores.index(max(scores))

    def _balance_partitions(self, partitions, graph, target_sizes):
        """
        Cân bằng kích thước các phân vùng bằng cách di chuyển các node
        """
        # Kiểm tra xem có cần cân bằng không
        need_balance = any(size < 0 for size in target_sizes)
        
        if not need_balance:
            return
            
        # Tìm phân vùng cần thêm node và phân vùng cần bớt node
        need_more = [i for i, size in enumerate(target_sizes) if size > 0]
        need_less = [i for i, size in enumerate(target_sizes) if size < 0]
        
        for from_idx in need_less:
            from_partition = partitions[from_idx]
            
            # Sắp xếp các node theo mức độ kết nối trong phân vùng (ít nhất lên đầu)
            nodes_to_move = sorted(from_partition, 
                                 key=lambda n: sum(1 for p_node in from_partition 
                                                 if graph.has_edge(n, p_node) or graph.has_edge(p_node, n)))
            
            # Số node cần di chuyển
            nodes_to_move = nodes_to_move[:abs(target_sizes[from_idx])]
            
            for node in nodes_to_move:
                # Tìm phân vùng tốt nhất để di chuyển node
                best_to = need_more[0]
                for to_idx in need_more:
                    if target_sizes[to_idx] > target_sizes[best_to]:
                        best_to = to_idx
                
                # Di chuyển node
                from_partition.remove(node)
                partitions[best_to].append(node)
                target_sizes[from_idx] += 1
                target_sizes[best_to] -= 1
                self.node_to_partition[node] = best_to
                
                # Cập nhật danh sách phân vùng cần thêm node
                if target_sizes[best_to] == 0:
                    need_more.remove(best_to)
                
                # Kiểm tra xem đã cân bằng xong chưa
                if target_sizes[from_idx] == 0:
                    break

    def identify_interfaces(self) -> List[Tuple[int, int, int, int]]:
        """
        Xác định các cạnh giao diện giữa các phân vùng
        """
        print("Identifying interfaces between partitions...")
        start_time = time.time()
        
        interfaces = []
        
        # Tìm các cạnh giữa các phân vùng
        for node in range(self.DFG.get_ndata()):
            if node not in self.node_to_partition:
                continue
                
            partition_id = self.node_to_partition[node]
            
            for operand_set in self.DFG.operands[node]:
                for operand in operand_set:
                    if operand < 0:  # Bỏ qua operand âm
                        continue
                    if operand not in self.node_to_partition:
                        continue
                        
                    operand_partition = self.node_to_partition[operand]
                    if operand_partition != partition_id:
                        interfaces.append((operand, node, operand_partition, partition_id))
        
        self.interfaces = interfaces
        interface_time = time.time() - start_time
        print(f"Identified {len(interfaces)} interfaces in {interface_time:.3f} seconds")
        
        return interfaces

    def create_sub_dfg(self, partition_idx):
        """
        Tạo một DFG con từ phân vùng và các giao diện liên quan
        """
        partition = self.partitions[partition_idx]
        
        # Tạo danh sách các node đầu vào của phân vùng từ giao diện
        input_nodes = set()
        for src, dst, src_part, dst_part in self.interfaces:
            if dst_part == partition_idx:
                input_nodes.add(src)
        
        # Tạo DFG mới
        sub_dfg = dfg.Dfg()
        
        # Xác định các node đầu vào
        for node in sorted(input_nodes):
            node_name = f"interface_{node}"
            sub_dfg.create_input(node_name)
        
        # Map từ node ID gốc sang node ID trong sub_dfg
        node_map = {}
        reverse_map = {}
        
        # Đánh dấu các node đầu vào
        for i, node in enumerate(sorted(input_nodes)):
            node_map[node] = i
            reverse_map[i] = node
        
        # Đánh dấu các node trong phân vùng
        for i, node in enumerate(sorted(partition)):
            node_map[node] = i + len(input_nodes)
            reverse_map[i + len(input_nodes)] = node
        
        # Tạo cấu trúc operands
        sub_dfg.ndata = len(input_nodes) + len(partition)
        sub_dfg.ninputs = len(input_nodes)
        sub_dfg.operands = [[] for _ in range(sub_dfg.ndata)]
        sub_dfg.exconds = [[] for _ in range(sub_dfg.ndata)]
        
        # Sao chép thông tin operands
        for node in partition:
            sub_node_idx = node_map[node] - sub_dfg.ninputs
            sub_dfg.oprtypes = self.DFG.oprtypes.copy() if hasattr(self.DFG, 'oprtypes') else []
            
            for operand_set in self.DFG.operands[node]:
                new_set = set()
                for operand in operand_set:
                    if operand < 0:
                        new_set.add(operand)  # Giữ nguyên các operand âm
                    elif operand in node_map:
                        new_set.add(node_map[operand])
                    else:
                        # Thêm node đầu vào ảo nếu không tìm thấy
                        virtual_input = sub_dfg.ninputs
                        sub_dfg.ninputs += 1
                        sub_dfg.ndata += 1
                        node_map[operand] = virtual_input
                        reverse_map[virtual_input] = operand
                        sub_dfg.create_input(f"virtual_{operand}")
                        new_set.add(virtual_input)
                
                sub_dfg.operands[sub_node_idx].append(new_set)
        
        # Đánh dấu các node đầu ra
        sub_dfg.output_ids = lambda: [node_map[n] - sub_dfg.ninputs for n in self.DFG.output_ids() if n in node_map and n in partition]
        
        return sub_dfg, node_map, reverse_map

    def map_partition(self, partition_idx, max_cycles=None):
        """
        Ánh xạ một phân vùng DFG lên CGRA
        """
        partition = self.partitions[partition_idx]
        if not max_cycles:
            max_cycles = min(len(partition) * 2, self.max_cycles // 2)
        
        print(f"\nMapping partition {partition_idx} with max_cycles={max_cycles}...")
        start_time = time.time()
        
        # Tạo sub-DFG
        sub_dfg, node_map, reverse_map = self.create_sub_dfg(partition_idx)
        
        # Lưu trữ trạng thái global từ main module
        original_DFG = main.DFG
        original_CGRA = main.CGRA
        original_nvars = main.nvars if hasattr(main, 'nvars') else 0
        original_all_clauses = main.all_clauses if hasattr(main, 'all_clauses') else []
        original_solver = main.solver if hasattr(main, 'solver') else None
        original_ncycles = main.ncycles if hasattr(main, 'ncycles') else 0
        
        # Thiết lập trạng thái mới cho main module
        main.DFG = sub_dfg
        main.CGRA = self.CGRA
        main.nvars = 0
        main.all_clauses = []
        main.solver = Glucose4(incr=True)
        
        # Thử ánh xạ với số chu kỳ tăng dần
        result = None
        for cycles in range(3, max_cycles + 1, 2):
            main.all_clauses = []
            main.nvars = 0
            main.ncycles = cycles
            
            print(f"  Trying with {cycles} cycles...")
            X = main.create_array(sub_dfg.ndata, self.CGRA.nnodes, cycles)
            Y = main.create_array(sub_dfg.ndata, len(self.CGRA.get_edges()), cycles)
            Z = main.create_array(sub_dfg.ndata, self.CGRA.nnodes, cycles)
            P = main.create_array(sub_dfg.nexs if hasattr(sub_dfg, 'nexs') else 0, self.CGRA.nnodes, cycles)
            Q = main.create_control_signal_array(sub_dfg)
            
            main.add_initial_conditions(X)
            main.add_final_conditions(X)
            main.add_existence_constraints(X, Y, Z)
            main.add_communication_constraints(X, Y)
            main.add_calculation_constraints(X, Y, Z, P, Q)
            main.add_block_constraints(P, Q, X, Y)
            main.add_control_signal_capacity(Q)
            main.add_capacity_constraints(X, Y, Z, self.nregs, self.nopes)
            
            if main.solver.solve():
                model = main.solver.get_model()
                result = {
                    'model': model,
                    'cycles': cycles,
                    'X': X,
                    'Y': Y,
                    'Z': Z,
                    'P': P,
                    'Q': Q,
                    'node_map': node_map,
                    'reverse_map': reverse_map,
                }
                print(f"  Success! Partition mapped in {cycles} cycles")
                break
            else:
                print(f"  Unsatisfiable with {cycles} cycles")
        
        # Khôi phục trạng thái global cho main module
        main.DFG = original_DFG
        main.CGRA = original_CGRA
        main.nvars = original_nvars
        main.all_clauses = original_all_clauses
        main.solver = original_solver
        main.ncycles = original_ncycles
        
        mapping_time = time.time() - start_time
        if result:
            print(f"Partition {partition_idx} mapped successfully in {result['cycles']} cycles ({mapping_time:.3f} seconds)")
        else:
            print(f"Failed to map partition {partition_idx} ({mapping_time:.3f} seconds)")
        
        return result

    def map_all_partitions(self):
        """
        Ánh xạ tất cả các phân vùng
        """
        self.local_mappings = []
        
        for i in range(len(self.partitions)):
            partition_size = len(self.partitions[i])
            max_cycles = min(partition_size * 3, self.max_cycles // 2)
            
            result = self.map_partition(i, max_cycles)
            
            if result:
                self.local_mappings.append(result)
            else:
                # Thử lại với max_cycles lớn hơn
                max_cycles *= 2
                result = self.map_partition(i, max_cycles)
                
                if result:
                    self.local_mappings.append(result)
                else:
                    print(f"Failed to map partition {i} even with {max_cycles} cycles")
                    return False
        
        return True

    def combine_mappings(self):
        """
        Kết hợp các ánh xạ phân vùng thành ánh xạ toàn cục
        """
        if not self.local_mappings:
            print("No local mappings to combine")
            return False
        
        print("\nCombining mappings...")
        start_time = time.time()
        
        # Xác định tổng số chu kỳ cần thiết
        total_cycles = max([result['cycles'] for result in self.local_mappings])
        
        # Thêm chu kỳ đệm để xử lý phụ thuộc giữa các phân vùng
        interface_overhead = len(self.interfaces) // 2
        total_cycles = min(total_cycles + interface_overhead, self.max_cycles)
        
        print(f"Attempting to combine with {total_cycles} cycles")
        
        # Khởi tạo SAT solver
        globals()['solver'] = Glucose4(incr=True)
        globals()['all_clauses'] = []
        globals()['nvars'] = 0
        globals()['ncycles'] = total_cycles
        
        # Tạo biến cho ánh xạ toàn cục
        X = create_array(self.DFG.ndata, self.CGRA.nnodes, total_cycles)
        Y = create_array(self.DFG.ndata, len(self.CGRA.get_edges()), total_cycles)
        Z = create_array(self.DFG.ndata, self.CGRA.nnodes, total_cycles)
        P = create_array(self.DFG.nexs, self.CGRA.nnodes, total_cycles)
        Q = create_control_signal_array(self.DFG)
        
        # Thêm các ràng buộc cơ bản
        add_initial_conditions(X)
        add_final_conditions(X)
        add_existence_constraints(X, Y, Z)
        add_communication_constraints(X, Y)
        add_calculation_constraints(X, Y, Z, P, Q)
        add_block_constraints(P, Q, X, Y)
        add_control_signal_capacity(Q)
        add_capacity_constraints(X, Y, Z, self.nregs, self.nopes)
        
        # Thêm ràng buộc từ ánh xạ cục bộ
        partition_to_cycle_offset = {}
        
        # Ánh xạ phân vùng vào các chu kỳ lần lượt
        current_cycle = 0
        for i, local_mapping in enumerate(self.local_mappings):
            partition_to_cycle_offset[i] = current_cycle
            current_cycle += local_mapping['cycles'] - 2  # Trừ đi chu kỳ đầu và cuối
        
        # Thêm các ràng buộc từ phân vùng
        for i, local_mapping in enumerate(self.local_mappings):
            cycle_offset = partition_to_cycle_offset[i]
            local_cycles = local_mapping['cycles']
            
            # Áp đặt vị trí của các node trong phân vùng
            for node in self.partitions[i]:
                local_node = local_mapping['node_map'][node]
                
                # Tìm vị trí của node trong local_mapping
                model_set = set(abs(x) for x in local_mapping['model'] if x > 0)
                
                for k in range(local_cycles):
                    for j in range(self.CGRA.nnodes):
                        if local_mapping['X'][local_node][j][k] in model_set:
                            # Map vào vị trí tương ứng trong global mapping
                            global_cycle = min(cycle_offset + k, total_cycles - 1)
                            add_clause([X[node][j][global_cycle]])
                            
                    for h in range(len(self.CGRA.get_edges())):
                        if local_mapping['Y'][local_node][h][k] in model_set:
                            global_cycle = min(cycle_offset + k, total_cycles - 1)
                            add_clause([Y[node][h][global_cycle]])
        
        # Thêm ràng buộc cho các cạnh giao diện
        for src, dst, src_part, dst_part in self.interfaces:
            # Src phải hoàn thành trước khi dst bắt đầu
            src_partition_offset = partition_to_cycle_offset[src_part]
            dst_partition_offset = partition_to_cycle_offset[dst_part]
            
            # Thêm ràng buộc truyền dữ liệu
            for k1 in range(total_cycles):
                for k2 in range(k1+1, total_cycles):
                    for j1 in range(self.CGRA.nnodes):
                        for j2 in range(self.CGRA.nnodes):
                            # Nếu src ở node j1 vào chu kỳ k1 và dst ở node j2 vào chu kỳ k2
                            # Đảm bảo dữ liệu có thể truyền từ j1 đến j2
                            # Đây là ràng buộc đơn giản hóa
                            add_clause([-X[src][j1][k1], -X[dst][j2][k2], Z[dst][j2][k2]])
        
        # Giải bài toán SAT kết hợp
        sat = globals()['solver'].solve()
        
        if sat:
            print("Global mapping found!")
            model = globals()['solver'].get_model()
            
            self.global_mapping = {
                'model': model,
                'X': X,
                'Y': Y,
                'Z': Z,
                'cycles': total_cycles
            }
            
            combining_time = time.time() - start_time
            print(f"Combining completed in {combining_time:.3f} seconds")
            return True
        else:
            print("Failed to find global mapping")
            
            # Thử lại với chu kỳ nhiều hơn nếu chưa đạt tối đa
            if total_cycles < self.max_cycles:
                new_total = min(total_cycles * 1.5, self.max_cycles)
                print(f"Retrying with {new_total} cycles")
                self.max_cycles = new_total
                return self.combine_mappings()
            
            combining_time = time.time() - start_time
            print(f"Combining failed after {combining_time:.3f} seconds")
            return False

    def solve(self):
        """
        Giải quyết vấn đề ánh xạ DFG lên CGRA sử dụng phương pháp phân cấp
        """
        start_time = time.time()
        
        # Bước 1: Phân vùng DFG
        self.partition_dfg()
        
        # Bước 2: Xác định giao diện
        self.identify_interfaces()
        
        # Bước 3: Ánh xạ từng phân vùng
        if not self.map_all_partitions():
            print("Failed to map all partitions")
            return False
        
        # Bước 4: Kết hợp ánh xạ
        if not self.combine_mappings():
            print("Failed to combine mappings")
            return False
        
        # Hiển thị kết quả
        model = self.global_mapping['model']
        X = self.global_mapping['X']
        Y = self.global_mapping['Y']
        Z = self.global_mapping['Z']
        
        # Interpret solution
        globals()['ncycles'] = self.global_mapping['cycles']
        interpret_solution(model, X, Y, Z)
        
        self.solve_time = time.time() - start_time
        print(f"Total hierarchical mapping time: {self.solve_time:.3f} seconds")
        return True

# Hàm wrapper để sử dụng từ main.py
def solve_mapping_hierarchical(num_cycles: int, nregs: int, nopes: int, num_partitions: int = 2):
    """
    Giải quyết vấn đề ánh xạ DFG lên CGRA sử dụng phương pháp phân cấp
    
    Args:
        num_cycles: Số chu kỳ tối đa
        nregs: Số thanh ghi cho mỗi PE
        nopes: Số operations mỗi PE có thể thực hiện trong một chu kỳ
        num_partitions: Số phân vùng cần chia
    """
    global DFG, CGRA
    
    # Kiểm tra DFG và CGRA đã được đặt từ bên ngoài
    if DFG is None or CGRA is None:
        if hasattr(main, 'DFG') and main.DFG is not None:
            DFG = main.DFG
        else:
            raise ValueError("DFG is None. Make sure to set it before calling solve_mapping_hierarchical")
        
        if hasattr(main, 'CGRA') and main.CGRA is not None:
            CGRA = main.CGRA
        else:
            raise ValueError("CGRA is None. Make sure to set it before calling solve_mapping_hierarchical")
    
    # Đảm bảo solver trong main module là một đối tượng Glucose4 hợp lệ
    if hasattr(main, 'solver') and main.solver is None:
        from pysat.solvers import Glucose4
        main.solver = Glucose4(incr=True)
    
    # Tiếp tục với thuật toán phân cấp
    mapper = HierarchicalMapper(DFG, CGRA, num_cycles, nregs, nopes, num_partitions)
    return mapper.solve()