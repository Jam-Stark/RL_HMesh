import os
import subprocess
import sys

# --- 配置 ---
mesh_dir = os.path.join("data", "mesh_files")
output_base_dir = "data"
executable_path = os.path.join(mesh_dir, "mesh2qhex.exe")

# 文件大小分类阈值 (字节)
SIZE_SMALL_MB = 1
SIZE_MEDIUM_MB = 10
SIZE_SMALL_BYTES = SIZE_SMALL_MB * 1024 * 1024
SIZE_MEDIUM_BYTES = SIZE_MEDIUM_MB * 1024 * 1024

# 输出目录名称
DIR_SMALL = os.path.join(output_base_dir, "qhex_small")
DIR_MEDIUM = os.path.join(output_base_dir, "qhex_medium")
DIR_LARGE = os.path.join(output_base_dir, "qhex_large")
# --- 配置结束 ---

def get_output_dir(file_size_bytes):
    """根据文件大小返回目标输出目录"""
    if file_size_bytes < SIZE_SMALL_BYTES:
        return DIR_SMALL
    elif file_size_bytes < SIZE_MEDIUM_BYTES:
        return DIR_MEDIUM
    else:
        return DIR_LARGE

def convert_files():
    """查找并转换 .mesh 文件"""
    if not os.path.isdir(mesh_dir):
        print(f"错误: 输入目录 '{mesh_dir}' 不存在。")
        sys.exit(1)

    if not os.path.isfile(executable_path):
        print(f"错误: 可执行文件 '{executable_path}' 不存在。")
        sys.exit(1)

    print(f"开始在 '{mesh_dir}' 目录中查找 .mesh 文件...")
    converted_count = 0
    skipped_count = 0
    error_count = 0

    for filename in os.listdir(mesh_dir):
        if filename.lower().endswith(".mesh"):
            input_filepath = os.path.join(mesh_dir, filename)
            
            try:
                file_size = os.path.getsize(input_filepath)
                target_dir = get_output_dir(file_size)
                
                # 创建目标目录（如果不存在）
                os.makedirs(target_dir, exist_ok=True)
                
                # 构建输出文件路径
                base_filename = os.path.splitext(filename)[0]
                output_filename = f"{base_filename}.Qhex"
                output_filepath = os.path.join(target_dir, output_filename)

                # 检查输出文件是否已存在
                if os.path.exists(output_filepath):
                    print(f"跳过: 输出文件 '{output_filepath}' 已存在。")
                    skipped_count += 1
                    continue

                print(f"转换: '{input_filepath}' ({file_size / 1024 / 1024:.2f} MB) -> '{output_filepath}'")
                
                # 构建并执行转换命令，并向其发送换行符以处理 "按任意键继续"
                command = [executable_path, input_filepath, output_filepath]
                # 使用 input='\n' 模拟按 Enter 键
                result = subprocess.run(command, capture_output=True, text=True, check=False, input='\n')

                # 检查返回码。注意：即使程序要求按键，它通常也会返回 0 表示成功
                if result.returncode != 0:
                    print(f"  警告/错误: 转换 '{filename}' 可能失败或有警告。返回码: {result.returncode}")
                    print(f"  标准输出: {result.stdout.strip()}")
                    print(f"  标准错误: {result.stderr.strip()}")
                    error_count += 1
                else:
                    # 检查文件是否被创建在了当前工作目录
                    expected_output_in_cwd = os.path.join(os.getcwd(), output_filename)
                    if os.path.exists(expected_output_in_cwd):
                        try:
                            # 将文件移动到目标目录
                            os.rename(expected_output_in_cwd, output_filepath)
                            print(f"  成功: '{output_filepath}' 已创建并移动。")
                            converted_count += 1
                        except OSError as e:
                            print(f"  错误: 无法将 '{expected_output_in_cwd}' 移动到 '{output_filepath}': {e}")
                            error_count += 1
                            # 尝试删除错误位置的文件
                            try:
                                os.remove(expected_output_in_cwd)
                            except OSError:
                                pass # 忽略删除错误
                    elif os.path.exists(output_filepath):
                         # 如果文件奇迹般地出现在了正确位置
                         print(f"  成功: '{output_filepath}' 已创建。")
                         converted_count += 1
                    else:
                        print(f"  错误: 转换命令执行成功，但在预期位置 '{output_filepath}' 或当前目录 '{expected_output_in_cwd}' 未找到输出文件。")
                        print(f"  标准输出: {result.stdout.strip()}")
                        print(f"  标准错误: {result.stderr.strip()}")
                        error_count += 1

            except FileNotFoundError:
                print(f"错误: 文件 '{input_filepath}' 未找到（可能在处理过程中被删除）。")
                error_count += 1
            except Exception as e:
                print(f"处理 '{filename}' 时发生未知错误: {e}")
                error_count += 1

    print("\n--- 转换完成 ---")
    print(f"成功转换: {converted_count}")
    print(f"跳过 (已存在): {skipped_count}")
    print(f"发生错误: {error_count}")

if __name__ == "__main__":
    convert_files()
