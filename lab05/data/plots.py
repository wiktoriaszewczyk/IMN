import matplotlib.pyplot as plt
import numpy as np

k = [16,8,4,2,1]
it = []
s = []

for i in range(5):
    
    fnameV = "V_k"+str(k[i])+".txt"
    dataV = np.genfromtxt(fnameV, delimiter='\t')
    x = dataV[:,0]
    y = dataV[:,1]
    V = dataV[:,2]

    fnameS = "S_k"+str(k[i])+".txt"
    dataS = np.genfromtxt(fnameS, delimiter='\t')
    it.append(dataS[:,0])
    s.append(dataS[:,1])

    # plot V

    # plt.figure()
    # plt.pcolor(x, y, V, cmap='hot', vmin=np.amin(V), vmax=np.amax(V), shading='auto')
    # plt.colorbar()
    # plt.title("V(x,y), k = "+str(k[i]))
    # plt.show()
    # print(np.transpose(dataV))

plt.figure()
for i in range(5):
    plt.plot(it[i], s[i], label="k = "+str(k[i]))

plt.title("S(it)")
plt.xlabel("it")
plt.ylabel("S(it)")
plt.grid()
plt.legend()
plt.show()


# plt.figure()
#     plt.tricontourf(xy_array[0],xy_array[1],Vn2, levels=np.linspace(min(Vn2),max(Vn2),999))
#     plt.colorbar(ticks=np.linspace(min(Vn2),max(Vn2),11))

#     plt.title("Releksacja globalna\n Zrelaksowany potencjał V(x,y), wG = "+str(wG))
#     plt.show()

# for i in range(len(S_array)):
#         plt.plot(it_array[i][1:],S_array[i][1:], label="wG = "+str(wG_array[i]))
#     plt.xscale("log")
#     plt.title("Metoda relaksacji globalnej\n Zmiana całki S(it)")
#     plt.xlabel("it")
#     plt.ylabel("S(it)")
#     plt.grid()
#     plt.legend()

#     plt.show()


