import numpy as np
alphas = np.empty(3)
betas  = np.empty(3)

alphas[0], betas[0] = 0.4940910158,	2.861432625

alphas[1], betas[1] =0.4742891967 , 	2.982852202

alphas[2], betas[2] = 0.4673632077 , 	2.992027469
print(alphas)
print(betas)    
print (alphas*betas)