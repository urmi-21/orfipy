import matplotlib.pyplot as plt
import sys
import seaborn 

#read data
with open(sys.argv[1]) as f:
    d=f.read().splitlines()
d.pop(0)

#datasets
tools=[]
#at=[]
#hs=[]
#bac=[]
#at_hpc=[]
#hs_hpc=[]
#bac_hpc=[]
getorf=[]
orfipy=[]
orfm=[]

for l in d:
    temp=l.split('\t')
    tools.append(temp[0])
    if temp[0]=='getorf':
        getorf=[float(x) for x in temp[1:]]
    elif temp[0]=='orfipy':
        orfipy=[float(x) for x in temp[1:]]
    else:
        orfm=[float(x) for x in temp[1:]]


print(tools,orfipy)


#plot
# Make the plot
# Set position of bar on X axis
barWidth = 0.25
r1 = [i for i in range(len(orfipy))]
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]

print(r1,r2,r2)

plt.bar(r2, getorf, color='#66c2a5', width=barWidth, edgecolor='white', label='getorf')
plt.bar(r3, orfm, color='#fc8d62', width=barWidth, edgecolor='white', label='orfm')
plt.bar(r1, orfipy, color='#8da0cb', width=barWidth, edgecolor='white', label='orfipy')

# Add xticks on the middle of the group bars
plt.xlabel('group', fontweight='bold')
#plt.xticks([r + barWidth for r in range(len(orfipy))], ['At', 'HS', 'Bac'])
# Create legend & Show graphic
plt.legend()
plt.show()
