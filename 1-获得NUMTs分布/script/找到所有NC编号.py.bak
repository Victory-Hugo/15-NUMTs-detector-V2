import re

def extract_nc_numbers(fasta_file):
    nc_numbers = {}
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                match = re.match(r'>((NC_\d+\.\d+))', line)
                if match:
                    nc_number = match.group(1)
                    chr_name_match = re.search(r'chromosome (\w+)', line, re.IGNORECASE)
                    if chr_name_match:
                        chr_name = chr_name_match.group(1)
                    else:
                        chr_name = nc_number  # 如果未找到染色体名称，则使用 NC 编号
                    if chr_name not in nc_numbers:
                        nc_numbers[chr_name] = []
                    nc_numbers[chr_name].append(nc_number)
    return nc_numbers

# 提供你的 GRCh38_latest_genomic.fna 文件的路径
fasta_file_path = '/mnt/d/REFERENCE_GRCH/GRCh37_latest_genomic.fna'

nc_numbers = extract_nc_numbers(fasta_file_path)

# 打印所有提取到的染色体名称和对应的 NC_ 编号
for chr_name, nc_list in nc_numbers.items():
    print(f"{chr_name}: {', '.join(nc_list)}")

