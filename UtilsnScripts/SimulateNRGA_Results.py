import numpy as np
import open3d as o3d
import copy
import os
import glob
import matplotlib.pyplot as plt
 
def parse_Transformations(trajectory_Registered):
    ####NOTE:### info --> <opcrit, MSE, RMSE, fps, totalIteration, totalDuration, M+N >
    lines = open(trajectory_Registered).readlines()
    matrices = {}   # <--- 4x4 transformation

    for i in range(0, int(len(lines)/4)):
       ll = lines[4*i+0 : 4*i+4]
       ss = ''.join(ll)
       mat = np.fromstring(ss, dtype=float, sep=' ').reshape(4,4)
       matrices[i] = mat
    return matrices

def parse_InfoTrans(trajectory_Registered):
    ####NOTE:### info --> <opcrit, MSE, RMSE, fps, totalIteration, totalDuration, M+N >
    lines = open(trajectory_Registered).readlines()
    tolFrames = int(len(open(trajectory_Registered).readlines()) / 5)
    matrices = {}   # <--- 4x4 transformation
    info = {}       # <--- (opcrit, mse, rmse, fps, totalIter, elapsedTime, #M+N points)
    for i in range(tolFrames):
       ll = lines[5*i+1 : 5*i+5]
       ss = ''.join(ll)
       mat = np.fromstring(ss, dtype=float, sep=' ').reshape(4,4)
       matrices[i] = mat
       el = np.fromstring(lines[5*i], dtype= float, sep = ' ').reshape(1, 7)
       info[i] = el
    return matrices, info

radius_normal = 1.0
iterFilesPath = "/media/ali/PortableSSD/Dropbox/02_NRGA/NRGA_Results/Supplementary_Material/Results_BunnyMissing1/qualityPlys/"
#sourceMesh = o3d.io.read_triangle_mesh("/home/ali/Dropbox/04_DATASETS/SUPER4PCS/kinect/stage/7.ply")
sourceMesh = o3d.io.read_point_cloud("/media/ali/PortableSSD/Dropbox/02_NRGA/NRGA_Results/Supplementary_Material/Results_BunnyMissing1/bunny_Y.ply")
targetPcd = o3d.io.read_point_cloud("/media/ali/PortableSSD/Dropbox/02_NRGA/NRGA_Results/Supplementary_Material/Results_BunnyMissing1/bunny_X.ply")
targetPcd.paint_uniform_color([1.0, 0.1, 0.1])
sourceMesh.paint_uniform_color([0.1, 0.5, 1.0])
#sourceMesh.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=radius_normal, max_nn=200))
targetPcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=radius_normal, max_nn=200))


srcClone = copy.deepcopy(sourceMesh)
trgtClone = copy.deepcopy(targetPcd)


class SimulateBySavedFiles:
    def __init__(self, src, trgt, iterFilesPath):
        self.frame_num = 0
        self.source = src
        self.target = trgt
        self.iterFilesPath = iterFilesPath
        self.saveImgsPath = os.path.join(self.iterFilesPath,  "frames")
        self.simulationFiles = os.listdir(iterFilesPath) #sorted(glob.glob(iterFilesPath + "*.ply"))

    def custom_draw_geometry_with_rotation(self, pcd1, pcd2):
        def rotate_view(vis):
            if self.frame_num >= len(self.simulationFiles)-1:
                frame_num = 0
                self.source = copy.deepcopy(srcClone)
                
            ctr = vis.get_view_control()
            ctr.rotate(0.0, 0.0)
            image = vis.capture_screen_float_buffer(False)            
            pcd1.points = o3d.io.read_point_cloud(os.path.join(self.iterFilesPath , "Template-" + str(self.frame_num+1) + ".ply")).points
            pcd1.paint_uniform_color([0.0, 0.0, 1.0])
            pcd1.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=radius_normal, max_nn=200))

            vis.update_geometry(pcd1)
            vis.poll_events()
            vis.update_renderer()
            if not os.path.exists(self.saveImgsPath):
                os.makedirs(self.saveImgsPath)
            plt.imsave(os.path.join(self.saveImgsPath, str(self.frame_num).zfill(4) + ".png"), np.asarray(image), dpi = 1)
            self.frame_num = self.frame_num + 1
            return False
            
        o3d.visualization.draw_geometries_with_animation_callback([pcd1, pcd2], rotate_view)


class SimulateBySavedTransformations:
    def __init__(self, src, trgt):
        self.frame_num = 1
        self.source = src
        self.target = trgt
    def custom_draw_geometry_with_rotation(self, pcd1, pcd2):
        def rotate_view(vis):
            if self.frame_num >= len(transHist)/4:
                frame_num = 1
                self.source = copy.deepcopy(srcClone)
                
            import matplotlib.pyplot as plt
            ctr = vis.get_view_control()
            ctr.rotate(0.0, 0.0)
            image = vis.capture_screen_float_buffer(False)
            
            #pcd1.transform(np.matmul(transHist[self.frame_num-1][:4, :4], np.linalg.inv(transHist[self.frame_num][:4, :4])))
            #pcd1.transform(np.linalg.inv(transHist[self.frame_num-1][:4, :4]))
            pcd1.transform(transHist[self.frame_num][:4, :4])
            plt.imsave("frames/"+ str(self.frame_num).zfill(4) + ".png",np.asarray(image), dpi = 1)
            self.frame_num = self.frame_num +1
            vis.update_geometry(pcd1)
            vis.poll_events()
            vis.update_renderer()
            return False
            
        o3d.visualization.draw_geometries_with_animation_callback([pcd1, pcd2], rotate_view)
        
        
o3d.visualization.draw_geometries([sourceMesh, targetPcd])
rotation_vec = np.ndarray((3,1),dtype=np.float64)
rotation_vec[0]=30
#pcd1.rotate(rotation_vec)
#pcd2.rotate(rotation_vec)
#pcd2.rotate(rotation_vec)

c = SimulateBySavedFiles(sourceMesh, targetPcd, iterFilesPath)
while True:
  c.custom_draw_geometry_with_rotation(c.source, c.target)
