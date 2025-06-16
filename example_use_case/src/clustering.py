
import pandas as pd
from sklearn.decomposition import PCA
import sys


input_file = sys.argv[1]
data= pd.read_csv(input_file)
x = data[["age", "weight", "blood_pressure", "average_heart_rate", "blood_sugar", "cholesterol"]]
pca = PCA(n_components=2)
X_pca = pca.fit_transform(x)


with open("explained_variance.txt", "w") as f:
  f.write(f"Explained variance ratio for each component: {pca.explained_variance_ratio_}\n")

dico = {}
x, y = X_pca[:, 0], X_pca[:, 1]
risk = list(data['risk'])
dico['x'] = x
dico['y'] = y
dico['risk'] = risk

df = pd.DataFrame.from_dict(dico)
df.to_csv('points.csv', index=False)

