import sys
import re
from collections import defaultdict
import statistics
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

def parse_trf_file(file_path):
    """
    Parse TRF annotation file

    Parameters:
        file_path (str): Path to the TRF annotation file

    Returns:
        dict: TRF annotation data organized by sample name
    """
    trf_data = {}
    current_sample = None

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Check if it's a sample identifier line
            if line.startswith('@'):
                current_sample = line[1:]  # Remove the leading '@'
                if current_sample not in trf_data:
                    trf_data[current_sample] = []
                continue

            # If no current sample, skip this line
            if current_sample is None:
                continue

            # Parse TRF annotation line
            parts = line.split()
            if len(parts) < 2:
                continue

            try:
                start = int(parts[0])
                end = int(parts[1])
                trf_data[current_sample].append((start, end))
            except ValueError:
                continue

    return trf_data

def check_coverage(trf_data, sample_name, query_start, query_end):
    """
    Check if any TRF annotation line can completely cover the given interval

    Parameters:
        trf_data (dict): Parsed TRF annotation data
        sample_name (str): Sample name
        query_start (int): Start position of the query interval
        query_end (int): End position of the query interval

    Returns:
        list: List of TRF annotation lines that can cover the query interval
    """
    if sample_name not in trf_data:
        print("No data")
        return []

    covering_intervals = []
    for start, end in trf_data[sample_name]:
        if start <= query_start and end >= query_end:
            covering_intervals.append((start, end))

    return covering_intervals

def calculate_coverage(trf_data, sample_name, target_start, target_end):
    if sample_name not in trf_data:
        print("No data")
        return []
    # Create a set to record all covered positions
    covered = set()
    target_length = target_end - target_start + 1
    
    # Iterate through all intervals
    for start, end in trf_data[sample_name]:
        # Calculate overlap with target interval
        overlap_start = max(start, target_start)
        overlap_end = min(end, target_end)
        
        # If there's overlap, add all covered positions
        if overlap_start <= overlap_end:
            covered.update(range(overlap_start, overlap_end + 1))
    
    # Calculate coverage percentage
    coverage_percentage = (len(covered) / target_length) * 100
    return coverage_percentage

# Count the main annotation type and MEI type
def get_main_annotation_type(mei_list):
    annotation_counts = {}
    mei_type_counts = {}

    for mei in mei_list:
        # Count annotation types
        annotation = mei['annotation'] or 'unknown'
        annotation_counts[annotation] = annotation_counts.get(annotation, 0) + 1

        # Count MEI types
        mei_type = mei['mei_type']
        mei_type_counts[mei_type] = mei_type_counts.get(mei_type, 0) + 1

    if annotation_counts and mei_type_counts:
        main_annotation = max(annotation_counts.items(), key=lambda x: x[1])
        main_mei_type = max(mei_type_counts.items(), key=lambda x: x[1])
        return main_annotation[0], main_annotation[1], annotation_counts, main_mei_type[0], main_mei_type[1], mei_type_counts
    return None, 0, {}, None, 0, {}

# Perform DBSCAN clustering analysis based on lengths and distances
def cluster_mei_arrays_dbscan(lengths, distances, eps=0.5, min_samples=3):
    """
    Perform DBSCAN clustering on MEI arrays based on length and distance
    eps: Neighborhood radius for determining point neighborhoods
    min_samples: Minimum number of samples in the neighborhood of a core point
    """
    if len(lengths) <= 2 or len(distances) <= 2:
        return None, None, None
    
    # Ensure both lists have the same length
    min_len = min(len(lengths), len(distances))
    lengths = lengths[:min_len]
    distances = distances[:min_len]
    
    # Prepare dataset
    data = np.column_stack((lengths, distances))
    
    # Standardize data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data)
    
    # Execute DBSCAN clustering
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    cluster_labels = dbscan.fit_predict(scaled_data)
    
    #print(cluster_labels)
    # Calculate center points for each cluster (excluding noise points with label -1)
    unique_labels = set(cluster_labels)
    if -1 in unique_labels:
        unique_labels.remove(-1)  # Remove noise point label
    
    centers = []
    for label in unique_labels:
        cluster_points = data[cluster_labels == label]
        center = np.mean(cluster_points, axis=0)
        centers.append(center)
    
    return cluster_labels, centers, scaler

def analyze_mei_periodicity(input_file, trf_data):
    """
    Analyze periodic occurrence patterns of MEI elements in repeat Masker output file
    """
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File not found {input_file}")
        return
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Skip the first three header lines
    data_lines = lines[3:]

    # Group by sample name
    sample_blocks = defaultdict(list)
    #print(sample_blocks)
    for line_num, line in enumerate(data_lines, start=4):
        line = line.strip()
        if not line:
            continue

        # Split columns, handling possible space separation
        columns = re.split(r'\s+', line)

        #print(len(columns))
        # Check if there are enough columns (at least 15)
        if len(columns) < 15:
            continue

        # Check if there are enough columns (at least 15)
        if len(columns) < 15:
            continue
        # Extract required fields
        try:
            sample_name = columns[4]  # 5th column
            start_pos = int(columns[5])  # 6th column
            end_pos = int(columns[6])   # 7th column
            annotation = columns[10]   # 11th column
        except (ValueError, IndexError):
            continue

        # Check if it's an MEI element
        mei_keywords = ['SINE', 'LINE', 'Retroposon', 'DNA', 'RC', 'LTR']
        is_mei = any(keyword in annotation for keyword in mei_keywords)

        if is_mei:
            mei_length = end_pos - start_pos + 1
            # Determine MEI type
            mei_type = 'unknown'
            for keyword in mei_keywords:
                if keyword in annotation:
                    mei_type = keyword
                    break

            #print(annotation,is_mei,start_pos,end_pos, mei_type, mei_length)
            sample_blocks[sample_name].append({
                'start': start_pos,
                'end': end_pos,
                'length': mei_length,
                'annotation': annotation,
                'mei_type': mei_type
            })

    # Analyze periodicity for each sample block
    periodic_samples = []

    # Create a list to store all data
    MEI_data_list = []

    for sample_name, mei_list in sample_blocks.items():
        if len(mei_list) < 3:  # Need at least 3 MEIs to detect periodicity
            continue
        
        #print(mei_list)
        # Sort by start position
        mei_list.sort(key=lambda x: x['start'])

        # Group by mei_type
        mei_by_type = defaultdict(list)
        for mei in mei_list:
            mei_by_type[mei['mei_type']].append(mei)

        #print(sample_name)
        # Perform periodicity analysis for each MEI type
        for mei_type, type_mei_list in mei_by_type.items():
            if len(type_mei_list) < 3:  # Each type needs at least 3 MEIs to detect periodicity
                continue
            #print("New classification:")
            #print(type_mei_list)
            # Calculate MEI lengths and distances
            lengths = [mei['length'] for mei in type_mei_list]
            positions = [mei['start'] for mei in type_mei_list]

            # Calculate distances (distance between adjacent MEIs)
            distances = []
            for i in range(1, len(positions)):
                distance = positions[i] - positions[i-1] - type_mei_list[i-1]['length']
                distances.append(distance)

            # Add average as the last element
            if distances:
                avg_distance = sum(distances) / len(distances)
                distances.append(avg_distance)
            
            # Ensure lengths and distances have the same length, as distances has one less element than lengths
            # If needed, adjust the lengths array to match the length of distances
            if len(lengths) > len(distances):
                lengths = lengths[:len(distances)]

            # Execute DBSCAN clustering analysis
            cluster_labels, centers, scaler = cluster_mei_arrays_dbscan(lengths, distances)
            type_mei_list_max = []
            if cluster_labels is not None:
                #print(f"DBSCAN clustering analysis result: cluster labels: {cluster_labels}")
                
                # Count the size of each cluster
                unique_labels, counts = np.unique(cluster_labels, return_counts=True)
                
                # Find the largest cluster (excluding noise points)
                non_noise_indices = [i for i, label in enumerate(unique_labels) if label != -1]
                if non_noise_indices:
                    max_cluster_idx = non_noise_indices[np.argmax(counts[non_noise_indices])]
                    max_cluster_label = unique_labels[max_cluster_idx]
                    max_cluster_size = counts[max_cluster_idx]
                    #print(f"\nLargest cluster is label {max_cluster_label}, containing {max_cluster_size} elements")
                    
                    # Create type_mei_list_max array, storing elements in the largest cluster
                    max_cluster_indices = np.where(cluster_labels == max_cluster_label)[0]
                    type_mei_list_max = [type_mei_list[i] for i in max_cluster_indices]
                    dis_mei_list_max = [distances[i] for i in max_cluster_indices]
                
            # If both length and distance are regular (coefficient of variation < 20%), consider it periodic
            if len(type_mei_list_max) >= 3:
                lengths = [mei['length'] for mei in type_mei_list_max]
                starts  = [mei['start'] for mei in type_mei_list_max]
                ends    = [mei['end'] for mei in type_mei_list_max]
                annotations = [mei['annotation'] for mei in type_mei_list_max]
                distances = dis_mei_list_max


                avg_length = statistics.mean(lengths)
                avg_distance = statistics.mean(distances)

                avg_distance_for_CV = max(avg_distance, 50)
                # Calculate coefficient of variation
                length_cv = statistics.stdev(lengths) / avg_length if avg_length > 0 else 0
                distance_stdev = statistics.stdev(distances)
                distance_cv = distance_stdev / avg_distance if avg_distance > 0 else 0

                if length_cv > 0.2 or (distance_stdev > 15 and  distance_cv > 0.2):
                    continue
                
                # Iterate through all data items and add to the list
                for i in range(len(lengths)):
                    # Create a dictionary for each data item
                    data_item = {
                        'sample_name': sample_name,
                        'length': lengths[i],
                        'annotation': annotations[i],
                        'start': starts[i],
                        'end': ends[i],
                        'distance': distances[i]
                    }
                    # Add dictionary to the list
                    MEI_data_list.append(data_item)

                cycles = len(lengths)  # Number of cycles
                main_annotation_type, main_count, all_counts, main_mei_type, main_mei_type_count, all_mei_type_counts = get_main_annotation_type(type_mei_list_max)
                
                # Check coverage interval
                # Calculate minimum start and maximum end
                start_pos = min([mei['start'] for mei in type_mei_list_max])
                end_pos = max([mei['end'] for mei in type_mei_list_max])
                # 
                # 'start': start_pos,
                
                covering_intervals = check_coverage(trf_data, sample_name, start_pos, end_pos)
                covering_persent = calculate_coverage(trf_data, sample_name, start_pos, end_pos)
                #print("covering_persent",covering_persent, "sample_name", sample_name)
                #print()

                periodic_samples.append({
                    'sample_name': sample_name,
                    'mei_type': mei_type,  # Add current analyzed MEI type
                    'cycles': cycles,
                    'avg_length': avg_length,
                    'avg_distance': avg_distance,
                    'main_annotation_type': main_annotation_type,
                    'main_annotation_count': main_count,
                    'main_mei_type': main_mei_type,  # New: main MEI type
                    'main_mei_type_count': main_mei_type_count,  # New: count of main MEI type
                    'length_cv': length_cv,      # New: length coefficient of variation
                    'distance_stdev': distance_stdev,   # New: distance standard deviation
                    'covering_intervals': covering_intervals,  # New: covering intervals
                    'start_pos': start_pos,
                    'end_pos': end_pos,
                    'covering_persent': covering_persent,
                })
        #break

    if periodic_samples:
        for sample in periodic_samples:
            print(f"Sample: {sample['sample_name']}", end=" | ")
            print(f"Main_MEI_type: {sample['main_annotation_type']}", end=" | ")
            print(f"Main_MEI_count: {sample['main_mei_type_count']}", end=" | ")
            print(f"Cycles: {sample['cycles']}", end=" | ")
            print(f"Length: {sample['avg_length']:.1f}bp", end=" | ")
            print(f"Distance: {sample['avg_distance']:.1f}bp", end=" | ")
            print(f"Length_CV: {sample['length_cv']:.3f}", end=" | ")
            print(f"Distance_SD: {sample['distance_stdev']:.3f}", end=" | ")
            print(f"Covering_intervals: {sample['covering_intervals']}", end=" | ")
            print(f"Start_position: {sample['start_pos']}", end=" | ")
            print(f"End_position: {sample['end_pos']}", end=" | ")
            print(f"Coverage: {sample['covering_persent']:.2f}%")
    # Output result sample_blocks
    print(f"{input_file} has {len(sample_blocks)} samples, of which {len(periodic_samples)} meet periodicity criteria")

    # Write data to file
    with open('mei_data.txt', 'w', encoding='utf-8') as f:
        for item in MEI_data_list:
            f.write(f"{item['sample_name']}\t{item['length']}\t{item['annotation']}\t{item['start']}\t{item['end']}\t{item['distance']}\n")

def process_file_list(list_file):
    try:
        with open(list_file, 'r') as f:
            files = [line.strip() for line in f if line.strip()]
            for file in files:
                print(f"\nProcessing file: {file}")
                analyze_mei_periodicity(file)
    except FileNotFoundError:
        print(f"Error: File list not found {list_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file list: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage:")
        print("  Single file: python mei_periodicity_analyzer.py -f <repeat_masker.out file>")
        print("  File list: python mei_periodicity_analyzer.py -l <file list>")
        sys.exit(1)

    option = sys.argv[1]
    input_file = sys.argv[2]
    trf_data = parse_trf_file(sys.argv[3])

    if option == "-f":
        #print(f"\nProcessing file: {input_file}")
        analyze_mei_periodicity(input_file, trf_data)
    elif option == "-l":
        process_file_list(input_file)
    else:
        print("Error: Invalid parameter option. Please use -f or -l")
        sys.exit(1)