
#Arciniega Arellano Andrés Eduardo
#Método de G.S.
from math import *

#************************************************************************************************************************
def ImprimeMatriz(dim,matriz):
    for i in range(0,dim):
        for j in range(0,dim):
            #imprime el numero redondead0 por 5 cifras
            print("|",round(matriz[i][j], 10), end=" |") 
        print()
#************************************************************************************************************************ 
#Método que realiza el producto de dos matrices
def ProductodeMatrices(dim,Matriz1, Matriz2, Matriz_Producto):
    for i in range(0,dim):
        for j in range(0,dim):
            for k in range(0,dim):
                Matriz_Producto[i][j] += Matriz1[i][k]*Matriz2[k][j]
#************************************************************************************************************************ 
#Método que obtiene la matriz transpuesta
def TranspuestadeMatriz(dim, Matriz_Trans, Matriz_Original):
    for i in range(0,dim):
        for j in range(0,dim):
            Matriz_Trans[i][j] = Matriz_Original[j][i]

#************************************************************************************************************************      
#Obtencion de Ki para el Metodo de G.S.
def Ki(dim, Vi, Qj):
    
    ki=0                    #variable que representa el valor buscado
    num=den=0 #Inicializamos el num y denominador de la division

    #obtiene el numerador
    for i in range (0, dim):
        num=num + (Vi[i]*Qj[i])

    #obtiene el denominador
    for i in range (0, dim):
        den=den + (Qj[i]*Qj[i])
    
    ki=num/den
    return ki
#************************************************************************************************************************
#Metodo de GramSchmidt
def GramSchmidt(dim, Matriz, Matriz_Ortogonal):
    #Para G.S es necesario ki,vi,qi,qj
    ki=0
    vi=[]
    qi=[]
    qj=[]
    #Se agregara otra variable que representa la operacion entre parentesis de G.S
    Operacion=[]

    #Se inicializa las diferentes listas con espacios vacios para evitar errores y considerar espacios vacios
    for i in range(0,dim):
        Operacion.append(0)
        vi.append(0)
        qi.append(0)
        qj.append(0)

    #Se inicializan dos contadores para manejar iteraciones dentro del método
    iter_qi=0
    iter_qj=1
    
    #Como primer paso se igualan la matriz ortogonal con la matriz dada
    for i in range(0,dim):
         Matriz_Ortogonal[i][0] = Matriz[i][0]
    #se incrementa qi
    iter_qi=1

    #Se deberá llenar la matriz ortogonal
    while(iter_qi<dim):  #while que se de no traspasar las dimensiones de la matriz

        #Obtencion de variables necesarias para el método
        
        for i in range (0, dim):        #Obtencion de vi
            vi[i] = Matriz[i][iter_qi]
        
        
        while(iter_qj<=iter_qi): #Mientras qj sea menor o igual a qi se obtiene la siguiente variable
            for i in range (0, dim):
                qj[i] = Matriz_Ortogonal[i][iter_qj-1]   #Obtencion de #qj
            iter_qj=iter_qj+1

            #Obtencion de: <vi,qj>/<qj,qj> = ki
            ki=Ki(dim,vi, qj)      
           
            for i in range (0, dim):
                Operacion[i] = Operacion[i] - (ki*qj[i])
        
        #Se hace la suma de Vi+Operacion
        for i in range (0, dim):
            qi[i] = vi[i]+Operacion[i]

        #Se agrega qi a la Matriz ortogonal
        for i in range(0,dim):
            Matriz_Ortogonal[i][iter_qi] = qi[i]
        
        #Para evitar problemas en el futuro se reestablecen todos los valores del método

        iter_qi=iter_qi+1
        iter_qj=1

        for i in range(0,dim):
            Operacion[i]=0
#************************************************************************************************************************
#Ortonormalizacion de la Matriz Ortogonal
def Ortonormalizacion(dim, Matriz_Ortogonal, Matriz_Ortonormal): 
    #Es necesaria la declaracion de algunas variables para el método
    ui=[]
    wi=[]
    Norma_wi=0
    #Y un contador
    cont=0

    #Se generan espacios en blanco en los vectores
    for i in range(0,dim):
        ui.append(0)
        wi.append(0)
    
    for i in range(0,dim):
        #Se calculan las variables para utilizar:  ui = wi / Norma_wi

        
        for b in range (0,dim):
            wi[b]=Matriz_Ortogonal[b][cont] #Obtencion de wi

        
        for b in range (0,dim):
            Norma_wi=Norma_wi+(wi[b]*wi[b]) #Obtencion de la Norma_wi
        Norma_wi=sqrt(Norma_wi)

        
        for b in range (0,dim):
            ui[b]=(wi[b])/Norma_wi  #Obtencion de ui

        #Se agrega ui a la matriz Ortonormal
        for b in range(0,dim):
            Matriz_Ortonormal[b][cont] = ui[b]

        #Se incrementa el contador y se limpia la norma para futuras iteraciones
        cont=cont+1
        Norma_wi=0
#************************************************************************************************************************
#Metodo que comprueba si es triangular superior (haya 0 por debajo de la diagonal)
#Si lo es regresa TRUE
#Si no lo es regresa FALSE
def TriangularSuperior( dim, Matriz):
    cont=0
    Es_o_Noes=True
    i=1
    while(i<dim):
        cont=cont+1
        for j in range(0,cont):
            if(Matriz[i][j]>0.000000000000001 or Matriz[i][j]<-0.000000000000001):
                Es_o_Noes=False
                break
        i=i+1
    return Es_o_Noes
#************************************************************************************************************************
#Implementacion de Algoritmo QR
def Algoritmo_QR(dim, Matriz,Matriz_Ortogonal,Matriz_Q,MatrizQT):

    #Primero se definen las variables necesarias para poder utilizar el algoritmo
    Matriz_triangular=False

    Matriz_R=[]
    Matriz_Ai_incrementada=[]
    Matriz_Ai=[]
    #Se inicia un contador en 0
    cont = 0

    #Se inicializa todo en 0
    for i in range(0,dim):
        Matriz_Ai.append([0]*dim)
        Matriz_Ai_incrementada.append([0]*dim)
        Matriz_R.append([0]*dim)

    #Primero se sabe que Ai = A
    for j in range(0,dim):
        for i in range(0,dim):
            Matriz_Ai[i][j] = Matriz[i][j]

    #Se repite el siguiente método hasta que sea una matriz triangular superior
    while(Matriz_triangular == False):
        
        #Se obtiene la Matriz_Q de Ai
        GramSchmidt( dim, Matriz_Ai, Matriz_Ortogonal)
        Ortonormalizacion( dim, Matriz_Ortogonal, Matriz_Q)

        #Se obtiene la Matriz_R de Ai
        #Se busca la matriz transpuesta de Q:  Q -> Qt
        TranspuestadeMatriz( dim, MatrizQT, Matriz_Q)

        #Se multiplica R=(Qt)(A)
        ProductodeMatrices( dim, MatrizQT, Matriz_Ai, Matriz_R)

        #Se muestran las iteraciones:
        print("\n A"+str(cont)+" = Q"+str(cont)+"x R "+str(cont))
        for i in range(0,dim):
            for j in range(0,dim):
                print("|",round(Matriz_Ai[i][j], 5), end="\t")
            print("\t = ",end='\t')
            for j in range(0,dim):
                print("|",round(Matriz_Q[i][j], 5), end=" | ") 
            print("",end='\t $')
            for j in range(0,dim):
                print("\t |",round(Matriz_R[i][j], 5), end=" |")
            print()

        #Se obtiene la matriz Ai+1 = (R)(Q)
        ProductodeMatrices(dim,Matriz_R, Matriz_Q, Matriz_Ai_incrementada)

        #Se comprueba si Matriz_Ai_incrementada es diagonal superior
        if(TriangularSuperior(dim, Matriz_Ai_incrementada)==True):
            #si es así salimos del ciclo
            Matriz_triangular=True
        #De lo contrario, se vuelve recursivo hasta que Ai+1 sea la nueva Ai    
        else:

            print("Iteracion " + str(cont) + "\n")
            cont = cont + 1
            
            for j in range(0,dim):
                for i in range(0,dim):
                    Matriz_Ai[i][j] = Matriz_Ai_incrementada[i][j]

                    #Se limpian las matrices para futuras iteraciones 
                    Matriz_R[i][j]=0
                    Matriz_Ai_incrementada[i][j]=0

        

    #Los eigenvalores son los elementos de la matriz triangular superior que se encuentran en la diagonal
    for i in range(0,dim):
        print(str(Matriz_Ai_incrementada[i][i]), ", ")

#************************************************************************************************************************
#**********************************************Inicio de Programa********************************************************

col = fil = 0
print ("Matriz n x n")
print("N representa el numero de Vectores")
n = int(input("ingrese cuantos vectores (n): "))
col = fil = n

#Creacion de todas las listas 
Matriz=[]
Matriz_Ortonormal =[]
Matriz_Ortogonal =[]
Matriz_Transpuesta=[]
Matriz_Resultado=[]

#Se generan las matrices en 0
for i in range(fil):
    Matriz.append([0]*col)
    Matriz_Ortogonal.append([0]*col)
    Matriz_Ortonormal.append([0]*col)
    Matriz_Transpuesta.append([0]*col)
    Matriz_Resultado.append([0]*col)
    
#Se rellena la Matriz con vectores acomodados por columna
for j in range(col):
    for i in range(fil):
        Matriz[i][j] = float(input("Ingresa la componente M_"+str(j+1)+str(i+1)+": "))

#Se imprime la Matriz Original
print("Matriz: ")
ImprimeMatriz(n,Matriz)

print ("\n Matriz_Ortogonal con GramSchmidt: ")
gs=GramSchmidt( n, Matriz, Matriz_Ortogonal)
ImprimeMatriz(n,Matriz_Ortogonal)

print("\n Nueva base Matriz_Ortonormal: ")
Ortonormalizacion( n, Matriz_Ortogonal,Matriz_Ortonormal)
ImprimeMatriz(n,Matriz_Ortonormal)

print("\n Eigenvalores: \n")
Algoritmo_QR(n,Matriz,Matriz_Ortogonal,Matriz_Ortonormal,Matriz_Transpuesta)

 