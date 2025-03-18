from utils.febfuncs import *
import subprocess
import numpy as np
import csv

# Constants
P = 106658  # Internal pressure in dyn/cm^2
E = 10000000  # Young's modulus in dyn/cm^2
nu = 0.49  # Poisson's ratio
R = 0.05
#t = 0.01
z_min = 0.45
z_max = 0.55

feb_nopress_file_path = "data/nopress_1mm.feb"
feb_press_file_path = "data/press_1mm.feb"
master_feb_press_file_path = "data/master_press.feb"
log_nopress_file_path = "data/nopress_1mm.log"
log_press_file_path = "data/press_1mm.log"
exe_path = "/dagon1/achitsaz/app/FEBioStudio/bin/febio4"

t_R_values = np.linspace(0.05, 0.3, 10)  # Range of thickness-to-radius ratios

# Output file
output_file = "stress_comparison_results_1mm.csv"

with open(output_file, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["t","R","t/R", "Weighted_Avg_Diff_Eigen1", "Weighted_Avg_Diff_Eigen2", "Weighted_Avg_Diff_Eigen3"])
    
    for t_R in t_R_values:
        t = t_R * R  # Compute thickness
        print(f'\n* Start to run for thickness : {t}\n')

        # Initinalized feb file and run FEBio for non-prestress condition
        modify_feb(feb_nopress_file_path, -P, t, E, nu)
        result = subprocess.run([exe_path, "-i", feb_nopress_file_path], capture_output=True, text=True)

        if result.stderr:
            print(f"Errors: {result.stderr}")

        # Extract FEBio stress field
        inflated_coordinate = extract_coordinate(log_nopress_file_path)
        nopre_stress_field = extract_stress_dict(extract_eigenvalue(log_nopress_file_path))


        # modify feb file and run FEBio for pre-stress condition on inflated geomentry
        modify_feb(feb_nopress_file_path, -P, t, E, nu,coordinates=inflated_coordinate,new_file_path=feb_press_file_path)
        replace_mesh_tag(feb_press_file_path, master_feb_press_file_path, output_feb_path=feb_press_file_path)
        modify_feb(feb_press_file_path, -P, t, E, nu)
        result = subprocess.run([exe_path, "-i", feb_press_file_path], capture_output=True, text=True)

        if result.stderr:
            print(f"Errors: {result.stderr}")

        # Extract FEBio stress field from pre-stress
        pre_stress_field = extract_stress_dict(extract_eigenvalue(log_press_file_path))
        centroid = extract_centroid(log_press_file_path)
        filtered_ids = [entry[0] for entry in centroid if z_min <= entry[3] <= z_max]

        # read grids data of inflated geometry
        nodes, elements = read_feb_file(feb_press_file_path)

        # Calculate and print the area 
        areas = []
        for elem in elements:
            node_ids = elem.node_ids  # Accessing the node IDs from the Element object
            if len(node_ids) == 4:  # Quad element
                p1, p2, p3, p4 = [nodes[node_id - 1] for node_id in node_ids]  # Adjust for 1-based index
                area = calculate_quad_area(p1, p2, p3, p4)
            elif len(node_ids) == 3:  # Triangular element
                p1, p2, p3 = [nodes[node_id - 1] for node_id in node_ids]  # Adjust for 1-based index
                area = calculate_triangle_area(p1, p2, p3)
            else:
                area = 0  # Handle cases where the element is neither tri3 nor quad4

            if area != 0:  # Only append non-zero areas
                areas.append(area)

        # Compare stress fields
        _, weighted_average = compare_stress_fields(nopre_stress_field, pre_stress_field, filtered_ids, areas)

        # Save results
        writer.writerow([t] + [R] + [t_R] + list(weighted_average))

print(f"Results saved to {output_file}")


