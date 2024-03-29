from ACANN import ACANN
from Database import Database
from torch.nn.modules.loss import KLDivLoss,L1Loss, SmoothL1Loss
from torch.optim import Adam,Rprop,Adamax, RMSprop,SGD,LBFGS
from torch.utils.data import DataLoader
import torch
import csv
import numpy as np
import math

from collections import deque


nb_data = 100000
first_epochs = int(5001)
last_epochs = int(10000)
last_loss = deque(maxlen=5000)
filename = 'checkpoint_archi_256_10batch.pth'
archi = [256,256]
print("Restarting ACANN for {} training data".format(nb_data))
model = ACANN(8,1024,archi,drop_p=0.09).double()
model.load_state_dict(torch.load(filename))
print("Model loaded")
# Import the data
train_data = Database(csv_target="../Data/A_training.csv",csv_input="../Data/G_training_reduced.csv",nb_data=nb_data).get_loader()
validation_data=Database(csv_target="../Data/A_validation.csv",csv_input="../Data/G_validation_reduced.csv",nb_data=int(nb_data*0.1)).get_loader()

trainloader = DataLoader(train_data,batch_size=int(nb_data/10),shuffle=True)
validationloader = DataLoader(validation_data,batch_size=int(nb_data))
print("Data Loaded")


# Define a function for computing the validation score
def validation_score(nn_model):
    nn_model.eval()
    val_error=L1Loss()
    with torch.no_grad():
        G_val,A_val=next(iter(validationloader))
        prediction=nn_model.forward(G_val)
        score=val_error(prediction,A_val)
    nn_model.train()
    return score.item()


#Define the loss
error = L1Loss()
#Define the optimizer
optimizer = Adam(model.parameters())
#RMSPRO 10 - 2e-3
#ADAM 10 - 1.2e-3

# Training parameters
step=-1
print_every = 500
print("Starting the training")
std=1e-9
nb_parameters=model.count_parameters()
# Training
for e in range(first_epochs,last_epochs):
    model.train()
    std /= 1.025
    #  Load a minibatch
    mean_loss=[]
    for G,A in trainloader:
        #G+=torch.randn(G.size()).normal_(mean=0,std=std).type(torch.cuda.DoubleTensor)
        step+=1
        # restart the optimizer
        optimizer.zero_grad()
        # compute the loss
        prediction = model.forward(G)
        loss = error(prediction,A)
        mean_loss.append(loss.item())
        # Compute the gradient and optimize
        loss.backward()
        optimizer.step()

        # Write the result
        if step % print_every == 0:
            last_loss.append(loss.item())
            print("Epoch {}/{} : ".format(e+1,last_epochs),
                  "Training MAE = {} -".format(loss.item()),
                  "Validation MAE = {}".format(validation_score(model)))
            torch.save(model.state_dict(),filename)
        
    with open("MAE_validation_"+str(nb_data)+"data_"+str(nb_parameters)+"parameters.csv",'a') as f:
        writer=csv.writer(f,delimiter=',')
        writer.writerow([str(e),str(validation_score(model))])
    with open("MAE_training_"+str(nb_data)+"data_"+str(nb_parameters)+"parameters.csv",'a') as f:
        writer=csv.writer(f,delimiter=',')
        writer.writerow([str(e),str(np.mean(mean_loss))])

print("Final score for {} data is : {}".format(nb_data,validation_score(model)))

