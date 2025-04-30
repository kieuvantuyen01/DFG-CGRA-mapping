"""
Tập tin cấu hình cho ứng dụng ánh xạ DFG lên CGRA
"""

# Tham số cơ bản
DEFAULT_NUM_CYCLES = 52  # Increased by 50% from 52
DEFAULT_NREGS = 2
DEFAULT_NOPES = 1

# Tham số cho phương pháp phân cấp
USE_HIERARCHICAL = True
NUM_PARTITIONS = 4  # Changed from 2 to 4 based on benchmark trend

# Tham số tối ưu hóa
FINCREMENT = True
FINTERPRET = True
FREDUCE = True
FMAXSAT = True  # Enable MaxSAT for better optimization

# Đường dẫn đến tập tin đầu vào
DFG_FILE = "f.txt"
CGRA_FILE = "e.txt"