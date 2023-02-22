# LOOP IN IMAGE

 Compile c++ code with CMake, or directly run the command line in the exe folder. 

* * *
```
  draw_surface_mesh.exe hh.hex 0 hh
``` 
Find the loop in hh.hex and save them to folder hh. The 0 here is only for debug, having no effect. 
***
```
draw_surface_mesh.exe leftbw1.hex 0 hh
```
Find the loop in leftbw1.hex (example in paper) and save them to folder hh.
***
**Note:** 

- This c++ code cannot read Dicom or nifti files directly, and you can use Matlab (untitled.mlx) to read them and convert them to the customized hex file.
- The visualization is only for debugging, not essentially part of our algorithm, and you can delete them in c++ code.

