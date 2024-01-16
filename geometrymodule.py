class geometry():
    def __init__(self, packagepath, trusspath):
        
        self.ppath = packagepath
        self.tpath = trusspath
        
        import os 
        import numpy as np
        import pandas as pd
        import math
        
        os.chdir(self.ppath)
        
        from utils import proputils as pu
        import main
        from names import GlobNames as gn
        
        os.chdir(self.tpath)
        
        n_nodes = np.linspace(0, 19, 20)
        xy_nodes = [[0, 0],
             [1, 0],
             [1, 1],
             [2, 0],
             [2, 1],
             [3, 0],
             [3, 1],
             [4, 0],
             [4, 1],
             [5, 0],
             [5, 1],
             [6, 0],
             [6, 1],
             [7, 0],
             [7, 1],
             [8, 0],
             [8, 1],
             [9, 0],
             [9, 1],
             [10, 0]]

        nodes_members = [[0, 1, 1, 0],
                 [0, 2, 1, 1],
                 [1, 2, 1, 2],
                 [2, 3, 1, 3],
                 [2, 4, 1, 4],
                 [3, 4, 1, 5],
                 [4, 5, 1, 6],
                 [4, 6, 1, 7],
                 [5, 6, 1, 8],
                 [6, 7, 1, 9],
                 [6, 8, 1, 10],
                 [7, 8, 1, 11],
                 [8, 9, 1, 12],
                 [9, 12, 1, 12],
                 [8, 10, 1, 13],
                 [9, 10, 1, 14],
                 [10, 12, 1, 13],
                 [11, 12, 1, 11],
                 [12, 14, 1, 10],
                 [11, 14, 1, 9],
                 [13, 14, 1, 8],
                 [14, 16, 1, 7],
                 [13, 16, 1, 6],
                 [15, 16, 1, 5],
                 [16, 18, 1, 4],
                 [15, 18, 1, 3],
                 [17, 18, 1, 2],
                 [18, 19, 1, 1],
                 [1, 3, 1, 0],
                 [3, 5, 1, 0],
                 [5, 7, 1, 0],
                 [7, 9, 1, 0],
                 [9, 11, 1, 0],
                 [11, 13, 1, 0],
                 [13, 15, 1, 0],
                 [15, 17, 1, 0],
                 [17, 19, 1, 0]]

        self.nodes = pd.DataFrame(data=[[float(i) for i in j] for j in xy_nodes], index=map(int, n_nodes), columns=["x", "y"])
        self.members = pd.DataFrame(data=nodes_members, columns=["Node 1", "Node 2", "no element", "crosssection"])
        self.props = pu.parse_file('bridge_frequency.pro')
        self.density = float(self.props["model"]["truss"]["density"])
        self.geomfile = "bridge_test.geom"

        self.props['init']['mesh']['file'] = self.geomfile
        self.FrequencyRequirements = [20, 40, 60] #Hz

    def rewrite_file(self):
        with open(self.geomfile, "w") as file:
            file.write('node: node #, x-coordinate; y-coordinate \n')
            file.write(self.nodes.to_csv(index=True, sep=' ', header=False))
            file.write('\n')
            file.write('member: Node #1; Node #2; number of elements; cross-section type \n')
            file.write(self.members.to_csv(index=False, sep=' ', header=False))

    def update_nodes(self, new_y, returns=False):

        if len(new_y) != 5:
            raise Exception("Length of new y vector not equal to 5")

        for i in range(5):
            self.nodes.loc[i * 2 + 2]["y"] = new_y[i]
            self.nodes.loc[20 - (i * 2 + 2)]["y"] = new_y[i]

        if returns == True:
            geometry.display_nodes(self)

    def update_areas(self, new_areas, returns=False):

        if len(new_areas) != 15:
            raise Exception("Length of new areas vector not equal to 15")

        self.props['model']['truss']['area'] = new_areas

        if returns == True:
            geometry.print_areas(self)

    def compute_mass(self, returns=False):
        import numpy as np
        
        n = len(self.members)
        dist = np.zeros(n)
        mass = np.zeros(n)

        if type(self.props["model"]["truss"]["area"]) == str:
            area = list(self.props["model"]["truss"]["area"][1:-1].split(","))
            area = [float(i) for i in area]
        else:
            area = self.props["model"]["truss"]["area"]

        for i in range(n):
            node1 = self.members.loc[i]["Node 1"]
            node2 = self.members.loc[i]["Node 2"]
            cs = int(self.members.loc[i]["crosssection"])

            x1 = self.nodes.loc[node1]["x"]
            y1 = self.nodes.loc[node1]["y"]

            x2 = self.nodes.loc[node2]["x"]
            y2 = self.nodes.loc[node2]["y"]

            dist[i] = np.round(np.sqrt((x1 - x2)**2 + (y1 - y2)**2), 2)
            mass[i] = np.round(area[cs] * self.density * dist[i], 2)

        self.dist = dist
        self.mass = mass.sum()

        if returns == True:
            geometry.print_mass(self)

    def run_FEM(self, check=False, returns=False):
        import os
        import numpy as np
        
        os.chdir(self.ppath)
        
        import main
        from names import GlobNames as gn
        
        os.chdir(self.tpath)
        
        geometry.rewrite_file(self)
        globdat = main.jive(self.props)
        self.EigenFrequencies = globdat[gn.EIGENFREQS][0:3]/2/np.pi

        geometry.print_frequencies(self)

        if check==True:
            geometry.check_eigenfrequencies(self)

        if returns==True:
            return self.EigenFrequencies

    def check_eigenfrequencies(self, returns=False):
        import numpy as np
        
        if np.sum(self.EigenFrequencies > self.FrequencyRequirements) == 3:
            self.requirements = True
        else:
            self.requirements = False
            
        if returns == True:
            return self.requirements

    def display_nodes(self):
        display(self.nodes)

    def display_members(self):
        display(self.members)

    def print_areas(self):
        print(self.props['model']['truss']['area'])

    def print_mass(self):
        print(f"the total mass of the structure is {self.mass:.3f} kg")

    def print_distance(self):
        print(self.dist)

    def print_frequencies(self):
        print(f'Smallest three natural frequencies: {self.EigenFrequencies} Hz')

    def return_frequencies(self):
        return self.EigenFreqs
        
    def return_variables(self):
        x1 = self.nodes.loc[2]["y"]
        x2 = self.nodes.loc[4]["y"]
        x3 = self.nodes.loc[6]["y"]
        x4 = self.nodes.loc[8]["y"]
        x5 = self.nodes.loc[10]["y"]
        
        x6 = self.props["model"]["truss"]["area"][0]
        x7 = self.props["model"]["truss"]["area"][1]
        x8 = self.props["model"]["truss"]["area"][2]
        x9 = self.props["model"]["truss"]["area"][3]
        x10 = self.props["model"]["truss"]["area"][4]
        x11 = self.props["model"]["truss"]["area"][5]
        x12 = self.props["model"]["truss"]["area"][6]
        x13 = self.props["model"]["truss"]["area"][7]
        x14 = self.props["model"]["truss"]["area"][8]
        x15 = self.props["model"]["truss"]["area"][9]
        x16 = self.props["model"]["truss"]["area"][10]
        x17 = self.props["model"]["truss"]["area"][11]
        x18 = self.props["model"]["truss"]["area"][12]
        x19 = self.props["model"]["truss"]["area"][13]
        x20 = self.props["model"]["truss"]["area"][14]
        
        return x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20
    
    def plot_structure(self):
        import matplotlib.pyplot as plt
        plt.figure(figsize=(12,4))

        for i in range(len(self.members)):
            node1 = self.members.loc[i]["Node 1"]
            node2 = self.members.loc[i]["Node 2"]

            x1 = self.nodes.loc[node1]["x"]
            y1 = self.nodes.loc[node1]["y"]

            x2 = self.nodes.loc[node2]["x"]
            y2 = self.nodes.loc[node2]["y"]
            plt.plot([x1, x2], [y1, y2], "b-")

        for i in range(len(self.nodes)):
            plt.plot(self.nodes.loc[i]["x"], self.nodes.loc[i]["y"], "ro")
            plt.annotate(i, (self.nodes.loc[i]["x"] + 0.1, self.nodes.loc[i]["y"]))

        plt.axis("scaled")
        plt.axis("off")
        plt.show()