import re
import shutil
import numpy as np
import xml.etree.ElementTree as ET
class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class Element:
    def __init__(self, node_ids):
        self.node_ids = node_ids

def read_feb_file(filename):
    # Parse the XML file
    tree = ET.parse(filename)
    root = tree.getroot()

    nodes = []
    elements = []

    # Find the Nodes section and parse node coordinates
    nodes_section = root.find(".//Nodes")
    for node in nodes_section.findall("node"):
        node_id = int(node.get('id'))
        coordinates = node.text.strip().split(',')
        x, y, z = float(coordinates[0]), float(coordinates[1]), float(coordinates[2])
        nodes.append(Point(x, y, z))

    # Find the Elements section and parse element connectivity
    elements_section = root.find(".//Elements")
    for elem in elements_section.findall("elem"):
        node_ids = list(map(int, elem.text.strip().split(',')))
        elements.append(Element(node_ids))

    return nodes, elements
def cross_product(v1, v2):
    return np.array([
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0]
    ])

def calculate_triangle_area(p1, p2, p3):
    # Create vectors from p1 to p2 and p1 to p3
    v1 = np.array([p2.x - p1.x, p2.y - p1.y, p2.z - p1.z])
    v2 = np.array([p3.x - p1.x, p3.y - p1.y, p3.z - p1.z])

    # Cross product of vectors v1 and v2
    cross_prod = cross_product(v1, v2)

    # Area is half the magnitude of the cross product
    return 0.5 * np.linalg.norm(cross_prod)

def calculate_quad_area(p1, p2, p3, p4):
    # Split the quadrilateral into two triangles: (p1, p2, p3) and (p1, p3, p4)
    area1 = calculate_triangle_area(p1, p2, p3)
    area2 = calculate_triangle_area(p1, p3, p4)

    # Total area is the sum of the two triangle areas
    return area1 + area2
def modify_feb(old_file_path, new_pressure, new_thickness, new_stiffness, new_poisson, coordinates=None, new_file_path=None):
    try:
        # Copy the original file to the new file path if provided
        if new_file_path:
            shutil.copy(old_file_path, new_file_path)
            print(f"File copied from {old_file_path} to {new_file_path}")
        else:
            new_file_path = old_file_path  # Overwrite the original file if new_file_path is not provided

        with open(new_file_path, 'r') as file:
            lines = file.readlines()
        
        new_lines = []
        inside_nodes_section = False
        nodes_header = None

        for line in lines:
            # Modify material properties
            line = re.sub(r'(<pressure\s*lc="1">)(.*?)(</pressure>)', r'\g<1>' + str(new_pressure) + r'\g<3>', line)
            line = re.sub(r'(<shell_thickness>)(.*?)(</shell_thickness>)', r'\g<1>' + str(new_thickness) + r'\g<3>', line)
            line = re.sub(r'(<E>)(.*?)(</E>)', r'\g<1>' + str(new_stiffness) + r'\g<3>', line)
            line = re.sub(r'(<v>)(.*?)(</v>)', r'\g<1>' + str(new_poisson) + r'\g<3>', line)
            if coordinates:
                # Detect <Nodes ...> section
                if re.match(r'<Nodes.*?>', line.strip()):  
                    inside_nodes_section = True
                    nodes_header = line.strip()
                    new_lines.append(nodes_header + "\n")  # Keep the original tag structure

                    # If new coordinates are provided, add them
                    if coordinates:
                        for point_id, x, y, z in coordinates:
                            new_lines.append(f'    <node id="{point_id}">{x},{y},{z}</node>\n')
                    continue  # Skip old nodes and wait for </Nodes>

                # Skip existing nodes inside <Nodes> section
                if inside_nodes_section:
                    if "</Nodes>" in line:
                        inside_nodes_section = False  # End of node section
                        new_lines.append(line)  # Add closing </Nodes> tag
                    continue  # Skip existing nodes

            new_lines.append(line)  # Keep all other lines unchanged

        # Write the modified lines to the new file
        with open(new_file_path, 'w') as file:
            file.writelines(new_lines)
        
        print(f"Successfully updated the file {new_file_path}")
        print(f"Updated pressure to {new_pressure}, thickness to {new_thickness}, stiffness to {new_stiffness}, and Poisson's ratio to {new_poisson}")
        if coordinates:
            print(f"Replaced node coordinates with the provided data.")

    except Exception as e:
        print(f"An error occurred: {e}")
def extract_stress_field(log_file):
    with open(log_file, 'r') as file:
        lines = file.readlines()

    time_found = False
    data_section = False
    stress_data = []

    for line in lines:
        # Check if we found "Time = 1.0"
        if re.search(r'Time\s*=\s*1', line):
            time_found = True
            continue  # Move to the next line
        
        # Look for "Data = element stresses" only if Time = 1.0 was found
        if time_found and re.search(r'Data\s*=\s*element stresses', line):
            data_section = True
            continue  # Move to the next line
        
        # Read stress data after finding "Data = element stresses"
        if data_section:
            parts = line.strip().split()
            if len(parts) == 7:  # Expecting 7 columns (ID, sx, sy, sz, sxy, syz, sxz)
                try:
                    element_id = int(parts[0])  # First column is element ID
                    sx, sy, sz, sxy, syz, sxz = list(map(float, parts[1:]))  # Convert stress values to float
                    stress_data.append((element_id, sx, sy, sz, sxy, syz, sxz))
                except ValueError:
                    continue  # Skip lines that can't be parsed
            else:
                # Break if we reach the end of the eigenvalue data section
                break
    return stress_data
def extract_coordinate(log_file):
    with open(log_file, 'r') as file:
        lines = file.readlines()

    time_found = False
    data_section = False
    coordinate_data = []

    for line in lines:
        # Check if we found "Time = 1.0"
        if re.search(r'Time\s*=\s*1', line):
            time_found = True
            continue  # Move to the next line
        
        # Look for "Data = coordinate" only if Time = 1.0 was found
        if time_found and re.search(r'Data\s*=\s*coordinate', line):
            data_section = True
            continue  # Move to the next line
        
        # Read coordinate of points after finding "Data = coordinate"
        if data_section:
            parts = line.strip().split()
            if len(parts) == 4:  # Expecting 4 columns (ID, x, y, z)
                try:
                    point_id = int(parts[0])  # First column is element ID
                    x, y, z = map(float, parts[1:])  # Convert coordinates to floats
                    coordinate_data.append((point_id, x, y, z))
                except ValueError:
                    continue  # Skip lines that can't be parsed
            else:
                # Break if we reach the end of the eigenvalue data section
                break
    return coordinate_data
def extract_centroid(log_file):
    with open(log_file, 'r') as file:
        lines = file.readlines()

    time_found = False
    data_section = False
    coordinate_data = []

    for line in lines:
        # Check if we found "Time = 1.0"
        if re.search(r'Time\s*=\s*1', line):
            time_found = True
            continue  # Move to the next line
        
        # Look for "Data = coordinate" only if Time = 1.0 was found
        if time_found and re.search(r'Data\s*=\s*element centroid position', line):
            data_section = True
            continue  # Move to the next line
        
        # Read coordinate of points after finding "Data = coordinate"
        if data_section:
            parts = line.strip().split()
            if len(parts) == 4:  # Expecting 4 columns (ID, x, y, z)
                try:
                    ele_id = int(parts[0])  # First column is element ID
                    x, y, z = map(float, parts[1:])  # Convert coordinates to floats
                    coordinate_data.append((ele_id, x, y, z))
                except ValueError:
                    continue  # Skip lines that can't be parsed
            else:
                # Break if we reach the end of the eigenvalue data section
                break

    return coordinate_data
def extract_eigenvalue(log_file):
    with open(log_file, 'r') as file:
        lines = file.readlines()

    time_found = False
    data_section = False
    eigenvalue_data = []
    for line in lines:
        # Check if we found "Time = 1.0"
        if re.search(r'Time\s*=\s*1', line):
            time_found = True
            continue  # Move to the next line
        
        # Look for "Data = coordinate" only if Time = 1.0 was found
        if time_found and re.search(r'Data\s*=\s*eigenvalue of Cauchy stress tensor', line):
            data_section = True
            continue  # Move to the next line
        
        # Read coordinate of points after finding "Data = coordinate"
        
        if data_section:
            parts = line.strip().split()
            if len(parts) == 4:  # Expecting 4 columns (ID, s1, s2, s3)
                try:
                    ele_id = int(parts[0])  # First column is element ID
                    s1, s2, s3 = map(float, parts[1:])  # Convert coordinates to floats
                    eigenvalue_data.append((ele_id, s1, s2, s3))
                except ValueError:
                    continue  # Skip lines that can't be parsed
            else:
                # Break if we reach the end of the eigenvalue data section
                break

    return eigenvalue_data
def replace_mesh_tag(source_feb_path, target_feb_path, output_feb_path=None):
    try:
        # Read source FEB file
        with open(source_feb_path, 'r') as file:
            source_content = file.read()

        # Extract the content within the <Mesh> tag
        source_mesh = re.search(r'<Mesh>(.*?)</Mesh>', source_content, re.DOTALL)
        if not source_mesh:
            raise ValueError('No <Mesh> tag found in source file.')

        source_mesh_content = source_mesh.group(1)

        # Read target FEB file
        with open(target_feb_path, 'r') as file:
            target_content = file.read()

        # Replace the <Mesh> tag in the target file
        modified_content = re.sub(r'<Mesh>.*?</Mesh>', f'<Mesh>{source_mesh_content}</Mesh>', target_content, flags=re.DOTALL)

        # Define the output path
        output_path = output_feb_path if output_feb_path else target_feb_path

        # Write the modified content to the output file
        with open(output_path, 'w') as file:
            file.write(modified_content)

        print(f'Mesh tag replaced successfully and saved to {output_path}')

    except Exception as e:
        print(f'Error: {e}')
def extract_stress_dict(stress_field):
    stress_dict = {entry[0]: entry[1:] for entry in stress_field}  # Assuming each entry is [element_id, sx, sy, sz, sxy, sxz, syz]
    return stress_dict
def compare_stress_fields(nopre_stress_field, pre_stress_field, element_ids, areas):
    comparison = {}
    weighted_differences = []
    total_area = 0  # To keep track of the total area for weighted average calculation

    print(f"Total number of element IDs to compare: {len(element_ids)}")

    for element_id, area in zip(element_ids, areas):
        nopre_stress = nopre_stress_field.get(element_id)
        pre_stress = pre_stress_field.get(element_id)

        if nopre_stress is not None and pre_stress is not None:
            nopre_stress = np.array(nopre_stress)
            pre_stress = np.array(pre_stress)

            # Calculate relative difference
            relative_difference = 100 * (pre_stress - nopre_stress) / (np.abs(nopre_stress) + 1e-8)
            comparison[element_id] = relative_difference

            # Append weighted difference and update total area
            weighted_differences.append(relative_difference * area)
            total_area += area
        else:
            print(f"Element {element_id} not found in one of the stress fields.")

    # Calculate average and standard deviation across all components, weighted by area
    if weighted_differences and total_area > 0:
        weighted_differences = np.array(weighted_differences)
        weighted_average = np.sum(weighted_differences, axis=0) / total_area
        weighted_std_deviation = np.sqrt(
            np.sum((weighted_differences - weighted_average)**2, axis=0) / total_area
        )

        print("\nWeighted Average Relative Difference:", weighted_average)
        print("Weighted Standard Deviation:", weighted_std_deviation)

    return comparison, weighted_average
