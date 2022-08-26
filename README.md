# scrapper

These python scripts are used for data analyzing in the SCRAPPY article

### Dependencies
Python libraries: FlowCytometryTools numpy pandas matplotlib seaborn

FlowCytometryTools doesn't work with python3.10, so use 3.9, 3.8 etc

### Draw heatmaps
1. Choose which heatmaps to draw and scales of heatmaps in the draw_parameters.txt
2. Choose which data to draw in the drugs_to_draw.txt
3. Launch draw_heatmaps.py 

### Add data
1. Insert names of drugs and paths to directories with .fcs files. Name and path should be separated with tab. Directory with .fcs files should have 3 sets of files: 
- 01-Well-{tube}.fcs for control (C0);
- 02-Well-{tube}.fcs for MIC/2; 
- 03-Well-{tube}.fcs for MIC
2. Launch add_data.py

Author: Fedor Ryabov, fedorrik1@gmail.com
