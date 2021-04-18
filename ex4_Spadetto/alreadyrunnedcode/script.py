import os
Ns = [50*i+1 for i in range(1,30)]
print(Ns)

optimizers=["", " -O1", " -O2", " -O3"] #vary optimizer
folders=["optn", "opt1", "opt2", "opt3"] #folder 4 optimizer
outfiles=["matmul.txt","matmul1.txt","matmul2.txt","matmul3.txt","matmul4.txt"] #files gen by ./test

 #check if the output dir is present
if not os.path.exists("./out"):
    os.makedirs("./out") #if not make it


for i in range(4): #vary optimization of compiler
    #select optimizer and folder
    cdir=folders[i]
    copt=optimizers[i]
    #compile
    os.system("gfortran"+copt+" test3.f90 -o test")
    if True :
        for n in Ns:
            #change dim of matrices
            print ("dim = " , n , ":")
            f = open("dim.txt", "w+")
            f.write(str(n))
            f.close()
            #run program
            os.system("./test")
    #once collected the data fit and plot everything
    os.system("gnuplot ./g_script.txt")

    #check if the output dir is present
    if not os.path.exists("./out/"+cdir):
        os.makedirs("./out/"+cdir) #if not make it

    #move everything into the respective folder
    for f in outfiles:
        os.system("mv ./out/"+f+" ./out/"+cdir )
    os.system("mv ./file.png ./out/"+cdir )
    os.system("mv ./fit.log ./out/"+cdir )
    #generate last plot
os.system("gnuplot ./g_script2.txt")
