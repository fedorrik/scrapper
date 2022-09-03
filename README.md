# scrapper

These python scripts are used for data analysis and visualization in the SCRAPPY article

### Dependencies
Python libraries: FlowCytometryTools numpy pandas matplotlib seaborn

FlowCytometryTools doesn't work with python3.10, so use 3.9, 3.8 etc

### Draw heatmaps
1. Choose which heatmaps to draw and the scales of heatmaps in the parameters.txt
2. Choose which data to draw in the drugs_to_draw.txt
3. Launch draw_heatmaps.py

### Add data
1. Choose PE-A threshold which will separate live and dead cells in the parameters.txt
2. Insert names of drugs and paths to directories with .fcs files. Name and path should be separated with tab. Directory with .fcs files should have 3 sets of files: 
- 01-Well-{tube}.fcs for control (C0);
- 02-Well-{tube}.fcs for MIC/2; 
- 03-Well-{tube}.fcs for MIC
3. Launch add_data.py
- If drug you trying to add already exist in the data file scrapper will ask you if you want to a.overwrite or b.skip it. Type a or b and press Enter to choose.

Author: Fedor Ryabov, fedorrik1@gmail.com
