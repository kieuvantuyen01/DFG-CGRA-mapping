import time
import copy
import os
import pandas as pd
import matplotlib.pyplot as plt
import dfg
import graph
import main
import sys
import traceback
from pysat.solvers import Glucose4

# Thêm đường dẫn hiện tại vào sys.path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

try:
    import hierarchical_mapping
    hierarchical_available = True
except ImportError:
    hierarchical_available = False
    print("Warning: Hierarchical mapping module not available. Install NetworkX to use this feature.")

def run_benchmark(dfg_files=["f.txt"], cgra_file="e.txt", num_cycles=52, nregs=2, nopes=1, num_partitions_range=[2, 3, 4]):
    """
    Chạy benchmark so sánh hiệu suất giữa các phương pháp ánh xạ
    """
    results = []
    
    for dfg_file in dfg_files:
        print(f"\n=== Benchmarking for DFG: {dfg_file} ===")
        
        # Kiểm tra tồn tại của file
        if not os.path.exists(dfg_file):
            print(f"Error: DFG file {dfg_file} does not exist")
            continue
        if not os.path.exists(cgra_file):
            print(f"Error: CGRA file {cgra_file} does not exist")
            continue
        
        # Đọc kiến trúc CGRA
        CGRA = graph.Graph()
        CGRA.create_node("mem", "_extmem")
        try:
            CGRA.read(cgra_file)
        except Exception as e:
            print(f"Error reading CGRA file: {e}")
            continue
        
        # Đọc DFG
        DFG = dfg.Dfg()
        try:
            DFG.read(dfg_file)
            MAC = True
            DFG.gen_operands(MAC, True)
        except Exception as e:
            print(f"Error reading or processing DFG file: {e}")
            continue
        
        # Xác minh rằng DFG và CGRA được khởi tạo đúng
        if not hasattr(DFG, 'ndata') or not hasattr(CGRA, 'nnodes'):
            print("Error: DFG or CGRA not properly initialized")
            continue
        
        # Chạy ánh xạ thông thường
        print("\n--- Regular Mapping ---")
        # Gán trực tiếp vào biến toàn cục trong main
        main.CGRA = CGRA
        main.DFG = DFG
        main.all_clauses = []
        main.nvars = 0
        main.ncycles = num_cycles
        
        # Quan trọng: Khởi tạo solver trước khi chạy solve_mapping
        # Phải đảm bảo main.py có biến toàn cục 'solver'
        if not hasattr(main, 'solver'):
            print("Warning: main module doesn't have 'solver' attribute. Adding it...")
            main.solver = Glucose4(incr=True)
        else:
            main.solver = Glucose4(incr=True)
            
        start_time = time.time()
        regular_sat = False
        try:
            regular_sat = main.solve_mapping(num_cycles, nregs, nopes)
        except Exception as e:
            print(f"Error in regular mapping: {e}")
            traceback.print_exc()
        regular_time = time.time() - start_time
        
        results.append({
            'dfg_file': dfg_file,
            'method': 'Regular',
            'partitions': 'N/A',
            'time': regular_time,
            'success': regular_sat is not False,
            'cycles': num_cycles
        })
        
        # Chạy ánh xạ phân cấp với số phân vùng khác nhau
        if hierarchical_available:
            for num_partitions in num_partitions_range:
                print(f"\n--- Hierarchical Mapping (partitions={num_partitions}) ---")
                
                # Cập nhật các biến toàn cục trong hierarchical_mapping
                hierarchical_mapping.DFG = DFG
                hierarchical_mapping.CGRA = CGRA
                
                # Đảm bảo trạng thái các biến toàn cục trong main là sạch
                main.nvars = 0
                main.all_clauses = []
                main.solver = Glucose4(incr=True)
                main.ncycles = num_cycles
                
                # Đảm bảo các biến cần thiết trong hierarchical_mapping
                if hasattr(hierarchical_mapping, 'nvars'):
                    hierarchical_mapping.nvars = 0
                if hasattr(hierarchical_mapping, 'all_clauses'):
                    hierarchical_mapping.all_clauses = []
                if hasattr(hierarchical_mapping, 'solver'):
                    hierarchical_mapping.solver = Glucose4(incr=True)
                
                start_time = time.time()
                hierarchical_sat = False
                try:
                    hierarchical_sat = hierarchical_mapping.solve_mapping_hierarchical(num_cycles, nregs, nopes, num_partitions)
                except Exception as e:
                    print(f"Error in hierarchical mapping: {e}")
                    traceback.print_exc()
                hierarchical_time = time.time() - start_time
                
                results.append({
                    'dfg_file': dfg_file,
                    'method': 'Hierarchical',
                    'partitions': num_partitions,
                    'time': hierarchical_time,
                    'success': hierarchical_sat is not False,
                    'cycles': num_cycles
                })
        
    # Tạo DataFrame với kết quả
    df = pd.DataFrame(results)
    
    # Hiển thị kết quả
    print("\n=== Benchmark Results ===")
    print(df)
    
    # Vẽ biểu đồ kết quả
    if not df.empty:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        for dfg_file in dfg_files:
            df_subset = df[df['dfg_file'] == dfg_file]
            
            # Regular mapping
            regular = df_subset[df_subset['method'] == 'Regular']
            if not regular.empty:
                ax.bar(f"{dfg_file} (Regular)", regular.iloc[0]['time'], color='blue')
            
            # Hierarchical mapping
            hierarchical = df_subset[df_subset['method'] == 'Hierarchical']
            for i, row in hierarchical.iterrows():
                ax.bar(f"{dfg_file} (Hierarchical p={row['partitions']})", row['time'], color='green')
        
        ax.set_ylabel('Time (seconds)')
        ax.set_title('Mapping Performance Comparison')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig('benchmark_results.png')
        plt.show()
        
        # Lưu kết quả vào CSV
        df.to_csv('benchmark_results.csv', index=False)
    
    return df

if __name__ == "__main__":
    dfg_files = ["f.txt"]  # Thêm nhiều DFG file nếu cần
    run_benchmark(dfg_files=dfg_files, num_partitions_range=[2, 3, 4])