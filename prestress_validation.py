from utils.febfuncs import *
import subprocess
import numpy as np
import csv

# Constants
#P = 106658  # Internal pressure in dyn/cm^2
E = 10000000  # Young's modulus in dyn/cm^2
nu = 0.49  # Poisson's ratio
R = 0.15
#t = 0.01
z_min = 1.45
z_max = 1.55

feb_nopress_file_path = "data/nopress_3mm.feb"
feb_press_file_path = "data/press_3mm.feb"
master_feb_press_file_path = "data/master_press.feb"
log_nopress_file_path = "data/nopress_3mm.log"
log_press_file_path = "data/press_3mm.log"
exe_path = "/dagon1/achitsaz/app/FEBioStudio/bin/febio4"
element_id_to_find = 3500  # Element ID to extract stress  for 3mm.feb
#element_id_to_find = 541  # Element ID to extract stress  for 2mm.feb
# element_id_to_find = 511  # Element ID to extract stress for 1mmfeb

t_R_values = np.linspace(0.05, 0.3, 10)  # Range of thickness-to-radius ratios
P_values = [79993,93325,106658,119989,133322,146654,159986]

# Output file
output_file = "stress_comparison_results_3mm.csv"

with open(output_file, mode="w", newline="") as file:
    writer = csv.writer(file)
    for P in P_values:

        P_mmhg = P * 0.00075006
        writer.writerow(["P",P_mmhg,"mmHg"])
        writer.writerow(["t","R","t/R", "Weighted_Avg_Diff_Eigen1", "Weighted_Avg_Diff_Eigen2", "Weighted_Avg_Diff_Eigen3","error_nopress_feb_vs_thin_t","error_press_feb_vs_thin_t"])
        
        for t_R in t_R_values:
            t = t_R * R  # Compute thickness
            print(f'\n* Start to run for thickness : {t}\n')

            # Compute analytical solutions for thin approximation
            stress_thin_t = (P * R) / t 

            # Initinalized feb file and run FEBio for non-prestress condition
            modify_feb(feb_nopress_file_path, -P, t, E, nu)
            result = subprocess.run([exe_path, "-i", feb_nopress_file_path], capture_output=True, text=True)

            if result.stderr:
                print(f"Errors: {result.stderr}")

            # Extract FEBio stress field
            inflated_coordinate = extract_coordinate(log_nopress_file_path)
            nopre_stress_field = extract_stress_dict(extract_eigenvalue(log_nopress_file_path))
            nopre_stress_tensor = extract_stress_field(log_nopress_file_path)
            
            # Get stress for the selected element
            stress_feb = None
            for data in nopre_stress_tensor:
                if data[0] == element_id_to_find:
                    stress_feb = list(data[1:])  # Extract stress components
                    print(stress_feb)
                    break

            if stress_feb is None:
                print(f"Warning: No stress data found for element {element_id_to_find}")
                continue

            # Extract Cartesian stress components from FEBio
            sx, sy, sz, sxy, syz, sxz = stress_feb  # FEBio gives stress in Cartesian

            # Convert to Cylindrical Stress Components According to the element_id_to_find element position
            stress_nopres_feb_t = sy  # Circumferential stress


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
            pre_stress_tensor = extract_stress_field(log_press_file_path)
            
            # Get stress for the selected element
            stress_feb = None
            for data in pre_stress_tensor:
                if data[0] == element_id_to_find:
                    stress_feb = list(data[1:])  # Extract stress components
                    print(stress_feb)
                    break

            if stress_feb is None:
                print(f"Warning: No stress data found for element {element_id_to_find}")
                continue

            # Extract Cartesian stress components from FEBio
            sx, sy, sz, sxy, syz, sxz = stress_feb  # FEBio gives stress in Cartesian

            # Convert to Cylindrical Stress Components According to the element_id_to_find element position
            stress_pres_feb_t = sy  # Circumferential stress

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
            # Calculate errors for Circumfrencial stress
            error_nopress_feb_vs_thin_t = abs((stress_nopres_feb_t - stress_thin_t) / stress_thin_t) * 100
            error_press_feb_vs_thin_t = abs((stress_pres_feb_t - stress_thin_t) / stress_thin_t) * 100

            # Save results
            writer.writerow([t, R, t_R] + list(weighted_average) + [error_nopress_feb_vs_thin_t,error_press_feb_vs_thin_t])

print(f"Results saved to {output_file}")


