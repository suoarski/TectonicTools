        
import bpy

bpy.ops.mesh.landscape_add(
		refresh=True,
		sphere_mesh=True,
		subdivision_x=512,
		subdivision_y=512,
		mesh_size=12,
		height=0.1,
		maximum=0.1,
		minimum=-0.1
		)
