import matplotlib.pyplot as plt
import numpy as np

print("hello")

import json

# Open and read the JSON file
with open('../output.json', 'r') as file:
    data = json.load(file)
data = data["data"]
# Print the data

x = []
y = []
for cor in data:
    x.append(cor[0])
    y.append(cor[1])




#print(dataarray)
#print(data["data"])
plt.plot(x,y)
plt.show()
