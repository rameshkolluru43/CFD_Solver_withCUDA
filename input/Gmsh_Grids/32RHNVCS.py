# Parse .geo file to find Physical Groups
def parse_geo_file(file_path):
    physical_groups = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if "Physical Curve" in line or "Physical Surface" in line:
                physical_groups.append(line.strip())
    return physical_groups

# Path to the .geo file
geo_file_path = "cylinder_mesh.geo"  # Replace with the path to your .geo file

# Parse and print physical groups
physical_groups = parse_geo_file(geo_file_path)
print("Physical Groups in .geo File:")
for group in physical_groups:
    print(group)