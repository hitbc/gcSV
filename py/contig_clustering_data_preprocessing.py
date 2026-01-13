
import re
import os
import argparse

def get_color_code(v_tag, colA, colB, colC, colD):
    """Determine color based on the given conditions"""
    is_non = "NON" in v_tag
    color_result = ''
    if not is_non:
        if colA != 0 and colB == 0:
            color_result = 'yellow'
        elif colB != 0:
            color_result = 'green'
        else:
            color_result = 'black'
    else:
        if colC != 0 and colD == 0:
            color_result = 'purple'
        elif colD != 0 and colC == 0:
            color_result = 'blue'
        elif colD != 0 and colC != 0:
            color_result = 'red'
        else:
            color_result = 'black'
    
    #print(color_result, is_non, v_tag, colA, colB, colC, colD)
    return color_result


def simple_cluster(numbers, threshold=20):
    """
    对数字列表进行简易聚类
    Args:
        numbers: 数字列表
        threshold: 聚类阈值，默认为50
    Returns:
        聚类后的列表，每个子列表代表一个聚类
    """
    if not numbers:
        return []
    
    # 对数字进行排序
    sorted_numbers = sorted(numbers)
    clusters = []
    current_cluster = [sorted_numbers[0]]
    
    for num in sorted_numbers[1:]:
        # 如果当前数字与上一个数字的距离在阈值内，加入当前聚类
        if num - current_cluster[-1] <= threshold:
            current_cluster.append(num)
        else:
            # 否则，保存当前聚类并开始新的聚类
            clusters.append(current_cluster)
            current_cluster = [num]
    
    # 添加最后一个聚类
    clusters.append(current_cluster)
    return clusters

def process_file(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    results = []
    current_sample = None
    data_lines = []

    for line in lines:
        line = line.strip()

        if line.startswith("Sample"):
            # If we have a previous sample, add it to results
            if current_sample is not None:
                # Add the sample header with count
                if len(data_lines) > 1:
                    data_lines = data_lines[1:]
                #results.append(f"Sample {current_sample} {len(data_lines)}")
                # Add all data lines for this sample
                results.extend(data_lines)
                data_lines = []

            # Extract sample info
            match = re.search(r'Sample (.+): File exists', line)
            if match:
                current_sample = match.group(1)
        elif line.startswith("AT_subregion"):
            # Parse data line
            # Extract ref info and sv_type_count
            ref_match = re.search(r'ref(.+?) :', line)
            sv_count_match = re.search(r'sv_type_count: \[(\d+), (\d+), (\d+), (\d+), (\d+)\]', line)
            cluster_match = re.search(r'cluster_size (\d+)', line)
            bp_match = re.search(r'type_VINV_GREEN_sv_base_absolute_BP_in_ref: \[(.*?)\]', line)
            cluster_ONLY_purple_match = re.search(r'cluster_size_with_ONLY_Insert_VNTR_svs: (\d+)', line)

            if bp_match:
                # 将匹配到的字符串转换为整数列表
                bp_list = [int(x) for x in bp_match.group(1).split(', ')] if bp_match.group(1) else []
            else:
                # 如果没有匹配到，设置为空数组
                bp_list = []

            if ref_match and sv_count_match and cluster_match and cluster_ONLY_purple_match:
                ref_info = ref_match.group(1)
                # Extract start, end, and V tag from ref_info
                # The format is chr1_101252500_101254500_0_2001_NON-VNTRFULL
                # We need to extract the start and end positions from the string
                # Use regex to extract the start and end positions
                #print("AA",ref_info)
                match = re.search(r'chr[0-9XYM]+_(\d+)_(\d+)_(\d+)+_(\d+)_([^_]+(?:_[^_]+)*)', ref_info)
                #refchr10_99837500_99840000_0_2501_NON-VNTRFULL
                #print("BB",match.group(1), match.group(2),match.group(3),match.group(4))
                
                if match:
                    start = match.group(3)
                    end = match.group(4)
                    v_tag = match.group(5)
                else:
                    # Fallback method
                    parts = ref_info.split('_')
                    # For format like chr1_101252500_101254500_0_2001_NON-VNTRFULL
                    # start is parts[-4], end is parts[-3], v_tag is parts[-2] + '_' + parts[-1]
                    if len(parts) >= 6:
                        start = parts[-4]
                        end = parts[-3]
                        v_tag = parts[-2] + '_' + parts[-1]
                    else:
                        start = parts[-4]
                        end = parts[-3]
                        v_tag = parts[-1]

                # Extract sv_type_count values
                colX, colA, colB, colC, colD = map(int, sv_count_match.groups())
                cluster_size = int(cluster_match.group(1))
                #print(cluster_ONLY_purple_match.group(1))
                cluster_size_only_purple = int(cluster_ONLY_purple_match.group(1))
                # Determine color
                color = get_color_code(v_tag, colA, colB, colC, colD)
                
                purple_red_cluster = 0
                if bp_list:
                    #print(bp_list)
                    clustered_list = simple_cluster(bp_list)
                    purple_red_cluster = len(clustered_list)
                    #print( current_sample ,"类别数", purple_red_cluster, "cluster_size_only_purple", cluster_size_only_purple, cluster_size, color, v_tag, "bp_list:",bp_list)  # 打印聚类结果
                    cluster_to_remove = cluster_size_only_purple - purple_red_cluster
                    cluster_size = cluster_size - max(cluster_to_remove, 0)
                    #print( current_sample ,"修正后： 类别数", purple_red_cluster, "cluster_size_only_purple", cluster_size_only_purple, cluster_size, color, v_tag, "bp_list:",bp_list)  # 打印聚类结果
                    

                # Add formatted data line
                data_lines.append(f"{current_sample} {start}_{end}_{v_tag} {color} {cluster_size} {colA} {colB} {colC} {colD}")

    # Add the last sample if exists
    if current_sample is not None:
        if len(data_lines) > 1:
            data_lines = data_lines[1:]
        #results.append(f"Sample {current_sample} {len(data_lines)}")
        results.extend(data_lines)

    # Write results to output file
    with open(output_file, 'w') as f:
        for result in results:
            f.write(result + '\n')

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process files for clustering analysis.')
    parser.add_argument('input_path', type=str, help='Path to the input file')
    parser.add_argument('--output', type=str, help='Path to the output file (optional)')
    args = parser.parse_args()
    
    input_path = args.input_path
    
    # If output path is not provided, use default naming
    if args.output:
        output_path = args.output
    else:
        # Extract file name from input path and add suffix
        file_name = os.path.basename(input_path)
        output_path = f"{file_name}.预处理1.txt"
    
    print(f"Processing {input_path} -> {output_path}")
    process_file(input_path, output_path)
    print(f"Completed processing {input_path}")

if __name__ == "__main__":
    main()
