from  torch.utils.data import TensorDataset
import pandas as pd
import torch
import numpy as np
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
class Database():
    def __init__(self, csv_target, csv_input, transform=None,nb_data=25000,nb_coefs=10):
        """
        Build the data set structure
        Args:
            csv_target (string): path to the target data (A(omega))
            csv_input (string) : path to the input data (G(tau))
            transform (callable, optional): Optional transform to be applied
                on a sample (eg: add noise).
        """
        self.input_data = pd.read_csv(csv_input,header=None,nrows=nb_data,usecols=[l for l in np.arange(nb_coefs)])
        self.target_data = pd.read_csv(csv_target,header=None,nrows=nb_data)
        self.transform = transform


    def get_loader(self):
        G=torch.tensor(self.input_data.values).double()
        A=torch.tensor(self.target_data.values).double()

        return TensorDataset(G.to(device),A.to(device))
