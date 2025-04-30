# Phương pháp phân cấp cho SAT-Based CGRA Mapping: Phương pháp luận và Chứng minh

## I. Cơ sở lý thuyết và Phương pháp luận

### A. Định nghĩa vấn đề
Vấn đề ánh xạ DFG lên CGRA có thể được định nghĩa chặt chẽ như sau:

**Cho:**
- DFG: G = (V, E) là đồ thị có hướng không chứa chu trình
- CGRA: C = (P, L) là kiến trúc gồm các PEs và đường kết nối
- Giới hạn tài nguyên: mỗi PE có nregs thanh ghi và nopes phép toán mỗi chu kỳ

**Tìm:**
- Ánh xạ M: V → P × T, với T là tập các chu kỳ
- Tối thiểu hóa |T| (số chu kỳ)

### B. Nguyên lý phân cấp (Hierarchical Approach)

Phương pháp phân cấp dựa trên ý tưởng chính:
1. **Phân vùng DFG thành các phần nhỏ hơn** 
2. **Giải quyết từng phần riêng biệt**
3. **Kết hợp lại các giải pháp cục bộ** thành giải pháp toàn cục

Cách tiếp cận này giải quyết vấn đề mở rộng quy mô của phương pháp SAT-based, vốn có độ phức tạp tăng theo hàm mũ của kích thước bài toán.

## II. Chứng minh tính đúng đắn

### A. Mô hình toán học 

Định nghĩa các biến quyết định:
- X(v,p,t) = 1 nếu node v được đặt tại PE p ở chu kỳ t
- Y(e,l,t) = 1 nếu cạnh e được ánh xạ lên đường kết nối l ở chu kỳ t
- Z(v,p,t) = 1 nếu node v được tính toán tại PE p ở chu kỳ t

### B. Tính đúng đắn của phương pháp phân cấp

**Định lý 1 (Tính đầy đủ):** *Nếu mỗi node trong DFG thuộc đúng một phân vùng và mỗi phân vùng được ánh xạ thành công, thì mọi node trong DFG đều được ánh xạ.*

**Chứng minh:**
- Phân vùng: V = V₁ ∪ V₂ ∪ ... ∪ Vₖ với Vᵢ ∩ Vⱼ = ∅, i≠j
- Mỗi v ∈ Vᵢ được ánh xạ thông qua giải pháp cục bộ của phân vùng i
- Vậy ∀v ∈ V, ∃p ∈ P, ∃t ∈ T: X(v,p,t) = 1

**Định lý 2 (Tính nhất quán):** *Nếu các ràng buộc phụ thuộc dữ liệu nội phân vùng được thỏa mãn và các ràng buộc liên phân vùng được đảm bảo trong giai đoạn kết hợp, thì mọi phụ thuộc dữ liệu trong DFG gốc đều được bảo toàn.*

**Chứng minh:**
- Với mọi cạnh (u,v) ∈ E trong cùng phân vùng, ánh xạ cục bộ đảm bảo t(v) > t(u)
- Với mọi cạnh (u,v) ∈ E giữa các phân vùng, giai đoạn kết hợp đảm bảo t(v) > t(u)
- Do đó, mọi phụ thuộc dữ liệu đều được thỏa mãn

**Định lý 3 (Độ phức tạp):** *Độ phức tạp thời gian của phương pháp phân cấp nhỏ hơn ánh xạ trực tiếp.*

**Chứng minh:**
- Độ phức tạp của SAT là O(2^n) với n là số biến
- Khi chia thành k phân vùng, mỗi phân vùng có khoảng |V|/k nodes
- Độ phức tạp của mỗi phân vùng là O(2^(n/k))
- Tổng độ phức tạp O(k·2^(n/k)) << O(2^n) khi k được chọn hợp lý

## III. Thuật toán phân cấp

### A. Thuật toán phân vùng

```
Function PartitionDFG(G, k):
    1. Tính toán mức (level) của mỗi node dựa trên đường đi từ đầu vào
    2. Gom các node theo mức thành các tập {L₁, L₂, ..., Lₘ}
    3. Khởi tạo k phân vùng rỗng G₁, G₂, ..., Gₖ
    4. Phân phối các node từ mỗi mức sao cho:
       a. Các node phụ thuộc được đặt trong cùng phân vùng khi có thể
       b. Cân bằng kích thước các phân vùng
       c. Tối thiểu số cạnh giữa các phân vùng
    5. Tối ưu hóa bằng thuật toán Kernighan-Lin
    6. Trả về {G₁, G₂, ..., Gₖ}
```

### B. Thuật toán xác định giao diện

```
Function IdentifyInterfaces(G, {G₁, G₂, ..., Gₖ}):
    1. Khởi tạo tập interfaces = ∅
    2. Với mỗi cạnh (u,v) ∈ E:
       a. Xác định phân vùng i chứa u và phân vùng j chứa v
       b. Nếu i ≠ j:
          i. Thêm (u,v,i,j) vào interfaces
    3. Trả về interfaces
```

### C. Thuật toán ánh xạ cục bộ

```
Function MapPartition(Gᵢ, C, max_cycles):
    1. Tạo sub-DFG với:
       a. Các node trong Vᵢ
       b. Các node đầu vào ảo cho các kết nối từ phân vùng khác
    2. Cho cycle = 1 đến max_cycles:
       a. Tạo biến SAT với cycle chu kỳ
       b. Thêm các ràng buộc SAT (tồn tại, giao tiếp, tính toán, tài nguyên)
       c. Giải bài toán SAT
       d. Nếu SAT, trả về ánh xạ (X_local, Y_local, Z_local, cycle)
    3. Trả về null (không tìm được ánh xạ)
```

### D. Thuật toán kết hợp ánh xạ

```
Function CombineMappings(interfaces, local_mappings, G, C, max_cycles):
    1. Khởi tạo các biến SAT cho ánh xạ toàn cục
    2. Thêm các ràng buộc từ các ánh xạ cục bộ:
       a. Các node đã được ánh xạ trong phân vùng phải giữ nguyên vị trí
    3. Thêm các ràng buộc đối với interfaces:
       a. Với mỗi cạnh (u,v) giữa các phân vùng, đảm bảo v được thực hiện sau u
    4. Thêm các ràng buộc tài nguyên toàn cục
    5. Giải bài toán SAT
    6. Trả về ánh xạ toàn cục (X, Y, Z) hoặc null
```

### E. Thuật toán phân cấp hoàn chỉnh

```
Function HierarchicalMapping(G, C, max_cycles, k):
    1. Phân vùng: partitions = PartitionDFG(G, k)
    2. Xác định giao diện: interfaces = IdentifyInterfaces(G, partitions)
    3. Ánh xạ từng phân vùng:
       a. local_mappings = []
       b. Với mỗi phân vùng partitions[i]:
          i. local_map = MapPartition(partitions[i], C, max_cycles/2)
          ii. Nếu local_map == null, tăng max_cycles/2 và thử lại
          iii. Thêm local_map vào local_mappings
    4. Kết hợp các ánh xạ:
       a. global_map = CombineMappings(interfaces, local_mappings, G, C, max_cycles)
       b. Nếu global_map == null, tăng max_cycles và thử lại
    5. Trả về global_map
```

## IV. Phân tích

### A. Các yếu tố ảnh hưởng đến hiệu suất

1. **Chất lượng phân vùng**: Phân vùng tốt với ít cạnh giữa các phân vùng dẫn đến ánh xạ hiệu quả
2. **Kích thước phân vùng**: Số lượng phân vùng k phải được tối ưu với kích thước DFG
3. **Chiến lược kết hợp**: Kết hợp hiệu quả các ánh xạ cục bộ là yếu tố quyết định thành công

### B. Phân tích độ phức tạp

- **Thời gian**: O(k·2^(n/k)) << O(2^n) với phương pháp thông thường
- **Không gian**: O(k·|V|·|P|·T_local) << O(|V|·|P|·T_global)

## V. Kết luận

Phương pháp phân cấp cho SAT-Based CGRA Mapping là một giải pháp đúng đắn về mặt lý thuyết để giải quyết vấn đề mở rộng quy mô. Phương pháp này vừa duy trì được tính chính xác của phương pháp SAT, vừa cho phép xử lý các DFG lớn hơn đáng kể.

Việc cài đặt phương pháp này đòi hỏi sự cẩn trọng, đặc biệt là trong các thuật toán phân vùng và kết hợp ánh xạ. Hiệu quả của phương pháp phụ thuộc mạnh mẽ vào các heuristic được sử dụng trong hai giai đoạn này.