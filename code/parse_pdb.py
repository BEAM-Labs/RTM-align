import pandas as pd 
import numpy as np
import argparse
import pdb

parser = argparse.ArgumentParser(description="parse pdb file")
parser.add_argument("--pdb_file_x", type=str)
parser.add_argument("--pdb_file_y", type=str)

def rmsd(coord1, coord2):
    n = len(coord2)
    d = coord1 - coord2
    return np.sqrt(np.sum(d * d) / n)

class BiPDBParser():
    def __init__(self, pdb_file_x, pdb_file_y, align_x, align_y, frag=8):
        # read initial pdb file
        self.df_x = pd.read_csv(pdb_file_x, header=None, sep="\s+")
        self.df_x.columns = ["ATOM", "atom_number", "type", "residual", "chain", "id", "x", "y", "z", "number1", "number2", "number3"]
        self.df_x[["x", "y", "z"]] = self.df_x[["x", "y", "z"]].astype(float)
        self.df_y = pd.read_csv(pdb_file_y, header=None, sep="\s+")
        self.df_y.columns = ["ATOM", "atom_number", "type", "residual", "chain", "id", "x", "y", "z", "number1", "number2", "number3"]
        self.df_y[["x", "y", "z"]] = self.df_y[["x", "y", "z"]].astype(float)
        # read alignment
        self.align_x = pd.read_csv(align_x, header=None, sep="\s+")
        self.align_x.columns = ["x", "y", "z"]
        self.align_x[["x", "y", "z"]] = self.align_x[["x", "y", "z"]].astype(float)
        self.align_y = pd.read_csv(align_y, header=None, sep="\s+")
        self.align_y.columns = ["x", "y", "z"]
        self.align_y[["x", "y", "z"]] = self.align_y[["x", "y", "z"]].astype(float)

        self.frag = frag
        
    def parse(self, save_name_x, save_name_y):
        # filter and map
        self.id_map_x, self.df_x = self.filter(self.df_x, self.align_x)
        self.id_map_y, self.df_y = self.filter(self.df_y, self.align_y)   

        # rotate and map
        self.rot_map = {}
        self.rotate()

        # rotate_x
        self.rotate_x()
        

        # self.df_x.dropna(inplace=True)
        # self.df_x[["chain", "id"]] = self.df_x[["chain", "id"]].astype(int)

        self.write(0, save_name_x)
        # self.write(1, save_name_y)

        
    
    def show(self, idx=0):
        if idx == 0:
            return self.df_x.head()
        else:
            return self.df_y.head()
    
    def write(self, idx, file_name):
        # rewrite df to pdb
        if idx == 0:
            atoms = [row for i, row in self.df_x.iterrows()]
        else:
            atoms = [row for i, row in self.df_y.iterrows()]
        formatted_atoms = "\n".join([
            "{:6} {:>4}  {:<4} {:>2} {:<1d}{:>4d} {:>11.3f} {:>7.3f} {:>7.3f}  {:<4.2f}  {:<4.2f} {:>11}"
            .format(atom["ATOM"], atom["atom_number"], atom["type"], atom["residual"], atom[4], atom[5], atom[6], atom[7], atom[8], atom[9], atom[10], atom[11])
            for atom in atoms
        ])
        # print(formatted_atoms)
        with open(file_name, "w") as f:
            f.writelines(formatted_atoms)
    
    def filter(self, pdbdf, align):
        # map alignment to original residuals
        prime = pdbdf[pdbdf["type"] == "C3'"]
        i = 0
        j = 0
        map = {}
        select_id = []
        while(j < len(align)):
            location = (align.iloc[j]["x"], align.iloc[j]["y"], align.iloc[j]["z"])
            while(i < len(prime)):
                location_origin = (prime.iloc[i]["x"], prime.iloc[i]["y"], prime.iloc[i]["z"])
                if location == location_origin:
                    map[j] = prime.iloc[i]["id"]
                    select_id.append(prime.iloc[i]["id"])
                    break
                i += 1
            j += 1
        pdbdf = pdbdf[pdbdf["id"].isin(select_id)]
        return map, pdbdf
    
    def rotate(self):
        align_to_rotate = {}
        # 与C++一样的截断代码
        i = 0
        Lali = len(self.align_x)
        while(i<Lali):
            if (Lali - i) <= int(self.frag * 1.5) :
                start_idx = i
                end_idx = Lali
                i = Lali
            elif i+self.frag > Lali :
                start_idx = i
                end_idx = Lali
                i = Lali
            else:
                start_idx = i
                end_idx = i + self.frag
                i = i + self.frag
            # 获得片段坐标
            frag_x = []
            frag_y = []
            for j in range(start_idx, end_idx):
                frag_x.append([self.align_x.iloc[j]["x"], self.align_x.iloc[j]["y"], self.align_x.iloc[j]["z"]])
                frag_y.append([self.align_y.iloc[j]["x"], self.align_y.iloc[j]["y"], self.align_y.iloc[j]["z"]])
            # 计算旋转矩阵
            rot_and_bias = self.kabsch(frag_x, frag_y)
            for j in range(start_idx, end_idx):
                align_to_rotate[j] = rot_and_bias
        self.rot_map = align_to_rotate

    def kabsch(self, movable_coords, reference_coords):
        movable_coords = np.array(movable_coords)
        reference_coords = np.array(reference_coords)
        c_movable        = np.mean(movable_coords, axis=0)
        c_reference      = np.mean(reference_coords, axis=0)
        # bias = c_movable - c_reference
        movable_coords   = movable_coords   - c_movable
        reference_coords = reference_coords - c_reference
        cm = np.dot(np.transpose(movable_coords), reference_coords)
        u, d, vt = np.linalg.svd(cm)
        rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
        
        return (rot, c_movable, c_reference)
    
    
    def rotate_x(self):
        id_to_rot = {}
        chain = 0
        last_rot = None
        for key, value in self.id_map_x.items():
            id_to_rot[value] = self.rot_map[key]
       
        for i in range(len(self.df_x)):
            id = self.df_x.iloc[i]["id"]
            rot, c_movable, c_reference = id_to_rot[id]
            location = np.array([self.df_x.iloc[i]["x"], self.df_x.iloc[i]["y"], self.df_x.iloc[i]["z"]]).reshape(1,3)
            location = location - c_movable 
            location = np.dot(location, rot)
            location = location + c_reference
            # if len(self.df_x) != 613:
            #     pdb.set_trace()
            
            self.df_x.iloc[i, 6] = location[0,0]
            self.df_x.iloc[i, 7] = location[0,1]
            self.df_x.iloc[i, 8] = location[0,2]
            if last_rot is None:
                last_rot = rot
            elif np.linalg.norm(last_rot - rot) > 0.0001:
                chain += 1
                last_rot = rot
            self.df_x.iloc[i, 4] = chain
            

args = parser.parse_args()
pdb1 = BiPDBParser(args.pdb_file_x, args.pdb_file_y, "output1.txt", "output2.txt")
pdb1.parse("RTMFAsup.pdb", "2et3_A_new.pdb")  
