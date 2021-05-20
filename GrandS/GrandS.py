import math     #importamos la libreria math para poder hacer uso de la funcion raiz cuadrada

matriz = []     #Declaramos nuestra matriz vacia
base_ortogonal = []
base_ortonormal = []
matriz_transpuesta = []
r = []
a0 = []
a1 = []
R = []

print ("Matriz mxm")
filas = int(input("Ingresa el valor de m: "))
#igualamos el numero de columnas con el numero de filas al ser una matrix mxm
columnas = filas

for i in range(filas):
    matriz.append([0]*columnas)
    base_ortogonal.append([0]*columnas)
    base_ortonormal.append([0]*columnas)
    matriz_transpuesta.append([0]*columnas)
    r.append([0]*columnas)
    a0.append([0]*columnas)
    a1.append([0]*columnas)
    R.append([0]*columnas)

# Pediremos los valores de la matriz por te clado
for j in range(columnas):
    for i in range(filas):
        matriz[i][j] = int(input("Ingresa la componente del vector "+str(j)+" posicion "+str(i)+": "))

# Funcion para obtener el vector
def obtenerVectores(matriz, filas, columnas):
    for i in range(columnas):
        print("Vector "+str(i)+":", end = " ")
        for j in range(filas):
            print(matriz[j][i], end = " ")
        print()

# Mandamos a llamar a la funcion para imprimir nuestros vectores
#obtenerVectores(matriz, filas, columnas)

# Funcion para imprimir la matriz
def imprimirMatriz(matriz, filas, columnas):
    for i in range(filas):
        for j in range(columnas):
            #redondeamos nuestros valores a 2 decimales
            print(round(matriz[i][j], 2), end=" ")
        print()

# Mandamos a llamar a la funcion para imprimir nuestra matriz
imprimirMatriz(matriz, filas, columnas)





def gramSchmidt(matriz, filas, columnas):
    # nuestra matriz va a ser igual a base_ortogonal para que podamos trabajar con ella
    for i in range(columnas):
        for j in range(filas):
            base_ortogonal[j][i] = matriz[j][i]
    
    for i in range(columnas):
        #recordemos que q0 = v0
        if i != 0: 
           productoPunto(matriz, i)

    for i in range(columnas):
        ortonormalizacion(base_ortogonal, i)
               
def productoPunto(matriz, i):
    cont=i  
    while(cont>0):
        prod_punt = div = division = escxvect = 0
        for j in range(filas):  
            prod_punt += matriz[j][i]*base_ortogonal[j][cont-1]
            div += base_ortogonal[j][cont-1]*base_ortogonal[j][cont-1]
        division = float(prod_punt/div)

        for j in range(filas):
            escxvect = division * base_ortogonal[j][cont-1]
            base_ortogonal[j][i] = base_ortogonal[j][i] - escxvect
        cont=cont-1

def ortonormalizacion(base_ortogonal, i):
    magnitud = 0
    for j in range(filas):
        magnitud += base_ortogonal[j][i]*base_ortogonal[j][i]
    raiz = math.sqrt(magnitud)
    
    for j in range(filas):
        base_ortonormal[j][i] = base_ortogonal[j][i]/raiz




def factorizacionQR(matriz, base_ortonormal):
    # A = QR  -> R = QtA
    # Primero sacaremos Qt
    matrizTranspuesta(base_ortonormal)
    # Por ultimo calcularemos R = Qt*A -> R = matriz_transpuesta*matriz
    multiplicacionMatrices(matriz_transpuesta, matriz)

def matrizTranspuesta(base_ortonormal):
    for i in range(filas):
        for j in range(columnas):
            matriz_transpuesta[i][j] = base_ortonormal[j][i]

def multiplicacionMatrices(matriz1, matriz2):
    for i in range(filas):
        for j in range(columnas):
            for k in range(filas):
                r[i][j] += matriz1[i][k]*matriz2[k][j]

def multiplicacionMatricesOrtogonales(matriz1, matriz2):
    for i in range(filas):
        for j in range(columnas):
            for k in range(filas):
                a1[i][j] += matriz1[i][k]*matriz2[k][j]




def algoritmoQR(a0, base_ortonormal):
    #Y DENTRO SERIA A0 = base_ortonormal*R  -> para sacar R
    factorizacionQR(a0, base_ortonormal)
    print("factorizacion QR")
    imprimirMatriz(r, filas, columnas)
    #Y LUEGO ES A1 = R*BASE_ORTONORMAL
    multiplicacionMatricesOrtogonales(r, base_ortonormal)
    print("Matriz A1")
    imprimirMatriz(a1, filas, columnas)

    #SI ES TRIANGULAR SUPERIOR A1, TERMINAMOS
    if(comprobacionTriangularSuperior(a1)):
        #imprimimos los valores de la diagonal principal
        for i in range(columnas):
            for j in range(filas):
                if(j == i):
                    print("Eigen valor "+str(i)+" es: "+str(a1[i][j]))
    else:
    #PERO SINO A0 = A1 Y SE REPITE EL ALGORITMO QR, GRAM_SCHMIDT Y FACTORIZACION QR
        for i in range(columnas):
            for j in range(filas):
                a0[j][i] = a1[j][i]

        gramSchmidt(a0, filas, columnas)
        algoritmoQR(a0, base_ortonormal)

def comprobacionTriangularSuperior(matrizA):
    bandera = True
    i = 1   #comenzamos desde la fila 1
    while (i<filas):
        j = 0
        while ((j<i) and (bandera == True)):
            if matrizA[i][j] != 0:
                bandera = False
                break
            else:
                j = j + 1
        i=i+1
    return bandera


gramSchmidt(matriz, filas, columnas)
print("Base ortogonal: ")
imprimirMatriz(base_ortogonal, filas, columnas)
print("Base ortonormal: ")
imprimirMatriz(base_ortonormal, filas, columnas)

#TAL VEZ AQUI PUEDA HACER A0 = MATRIZ
for i in range(columnas):
    for j in range(filas):
        a0[j][i] = matriz[j][i]

#imprimirMatriz(a0, filas, columnas)
algoritmoQR(a0, base_ortonormal)

#factorizacionQR(matriz, base_ortonormal)
#print("Matriz Transpuesta: ")
#imprimirMatriz(matriz_transpuesta, filas, columnas)
#print("Matriz R")
#imprimirMatriz(r, filas, columnas)

