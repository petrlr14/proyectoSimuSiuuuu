//check
float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}
/*CALCULOS PARA CHI*/
//check
void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);
}
//check
void calculateBeta(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) =  1; 
    B.at(0).at(6) =  0; 
    B.at(0).at(7) =  0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) =  1; 
    B.at(0).at(10) =  0; 
    B.at(0).at(11) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0; 
    B.at(1).at(6) = 1;
    B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1;
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
    B.at(2).at(4) = -1;
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0;
    B.at(2).at(7) = 1;
    B.at(2).at(8) = -1;
    B.at(2).at(9) = 0;
    B.at(2).at(10) = 0;
    B.at(2).at(11) = 1;
}
//check
void calculateGamma(Matrix &m, mesh msh, int eIndex){
	zeroes(m,12,3);
    element e = msh.getElement(eIndex);

    node node1=msh.getNode(e.getNode1()-1);
    node node2=msh.getNode(e.getNode2()-1);
    node node3=msh.getNode(e.getNode3()-1);
    node node4=msh.getNode(e.getNode4()-1);
    
    float x1, x2, x3, x4,y1, y2, y3, y4;
    
    x1= node1.getX();
    x2= node2.getX();
    x3= node3.getX();
    x4= node4.getX();

    y1= node1.getY();
    y2= node2.getY();
    y3= node3.getY();
    y4= node4.getY();

    float first, second, third, fourth;
    first =(float) (18*x1+9*x2+9*x3+9*x4+(2*(3*pow(y1,2)+(2*y1*(y2+y3+y4))+pow(y2,2)+(y2*(y3+y4))+pow(y3,2)+(y3*y4) +pow(y4, 2) )) )/360;
    
    second = (float)(9*x1+18*x2+9*x3+9*x4+(2*(pow(y1,2)+(y1*(2*y2+y3+y4))+3*pow(y2,2)+(y2*(y3+y4))+pow(y3,2)+(y3*y4) +pow(y4, 2) )) )/360;

    third =(float) (9*x1+9*x2+18*x3+9*x4+(2*(pow(y1,2)+(y1*(y2+2*y3+y4))+pow(y2,2)+(y2*(2*y3+y4))+3*pow(y3,2)+2*(y3*y4) +pow(y4, 2) )) )/360;

    fourth = (float)(9*x1+9*x2+9*x3+(2*(9*x4+pow(y1,2)+(y1*(y2+y3+2*y4))+pow(y2,2)+(y2*(y3+y4))+pow(y3,2)+2*(y3*y4) +3*pow(y4, 2) )) )/360;

	m.at(0).at(0) = first;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = second;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = third;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = fourth;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = first;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = second;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = third;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = fourth;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = first;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = second; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = third;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = fourth; 
	
}

/*CALCULOS PARA ZETA*/
//check
float calculateE(mesh msh, int eIndex){
    element e = msh.getElement(eIndex);

    node node1=msh.getNode(e.getNode1()-1);
    node node2=msh.getNode(e.getNode2()-1);
    node node3=msh.getNode(e.getNode3()-1);
    node node4=msh.getNode(e.getNode4()-1);
    
    float x1, x2, x3, x4,y1, y2, y3, y4;
    
    x1= node1.getX();
    x2= node2.getX();
    x3= node3.getX();
    x4= node4.getX();

    y1= node1.getY();
    y2= node2.getY();
    y3= node3.getY();
    y4= node4.getY();

    return (float)(pow(x1, 2)*(3*y1+y2+y3+y4)+x1*(pow(x2, 2)*(2*y1+2*y2+y3+y4)+x3*(2*y1+y2+y3+y4)+x4*(3*y1+y2+y3+y4))+pow(x2, 2)*(y1+3*y2+y3+y4)+x2*(x3*(y1+2*y2+2*y3+y4) +x4*(y1+2*y2+y3+2*y4)) + pow(x3, 2)*(y1+y2+3*y3+y4)+x3*x4*(y1+y2+2*(y3+y4)) + pow(x4,2)*(y1+y2+y3+3*y4)+160 )/360;
}


/*CALCULOS PARA THETA*/
//check
void calculateA(Matrix &m, mesh msh, int eIndex){
    zeroes(m, 3, 4);
    m.at(0).at(0) = -1;
    m.at(0).at(1) = 1;
    m.at(0).at(2) = 0;
    m.at(0).at(3) = 0;

    m.at(1).at(0) = -1;
    m.at(1).at(1) = 0;
    m.at(1).at(2) = 1;
    m.at(1).at(3) = 0;

    m.at(2).at(0) = -1;
    m.at(2).at(1) = 0;
    m.at(2).at(2) = 0;
    m.at(2).at(3) = 1;
}
//check
void calculateQ(Matrix &m, mesh msh, int eIndex){
    zeroes(m, 12, 3);
    element e = msh.getElement(eIndex);

    node node1=msh.getNode(e.getNode1()-1);
    node node2=msh.getNode(e.getNode2()-1);
    node node3=msh.getNode(e.getNode3()-1);
    node node4=msh.getNode(e.getNode4()-1);
    
    float y1, y2, y3, y4;
    
    y1= node1.getY();
    y2= node2.getY();
    y3= node3.getY();
    y4= node4.getY();

    float first, second, third, fourth;

    first = (float)7*(2*y1+y2+y3+y4)/120;
    second = (float)7*(y1+2*y2+y3+y4)/120;
    third = (float)7*(y1+y2+2*y3+y4)/120;
    fourth = (float)7*(y1+y2+y3+2*y4)/120;

    m.at(0).at(0) = first;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = second;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = third;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = fourth;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = first;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = second;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = third;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = fourth;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = first;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = second; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = third;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = fourth; 

}

/*CALCULOS PARA Pi*/
//check
void calculateH(Matrix &m, mesh msh, int eIndex){
    zeroes(m, 12, 3);
    element e = msh.getElement(eIndex);

    node node1=msh.getNode(e.getNode1()-1);
    node node2=msh.getNode(e.getNode2()-1);
    node node3=msh.getNode(e.getNode3()-1);
    node node4=msh.getNode(e.getNode4()-1);
    
    float z1, z2, z3, z4;
    
    z1= node1.getZ();
    z2= node2.getZ();
    z3= node3.getZ();
    z4= node4.getZ();

    float first, second, third, fourth;

    first = (float)(2*z1+z2+z3+z4)/40;
    second = (float)(z1+2*z2+z3+z4)/40;
    third = (float)(z1+z2+2*z3+z4)/40;
    fourth = (float)(z1+z2+z3+2*z4)/40;

    m.at(0).at(0) = first;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = second;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = third;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = fourth;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = first;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = second;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = third;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = fourth;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = first;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = second; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = third;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = fourth; 

}

/*Calculos para Tao*/
void calculateO(Matrix &m, mesh msh, int eIndex){
    zeroes(m, 4, 3);
    element e = msh.getElement(eIndex);

    node node1=msh.getNode(e.getNode1()-1);
    node node2=msh.getNode(e.getNode2()-1);
    node node3=msh.getNode(e.getNode3()-1);
    node node4=msh.getNode(e.getNode4()-1);
    
    float z1, z2, z3, z4;
    
    z1= node1.getZ();
    z2= node2.getZ();
    z3= node3.getZ();
    z4= node4.getZ();

    float first, second, third, fourth;

    first = (float)(2*z1+z2+z3+z4)/120;
    second =(float) (z1+2*z2+z3+z4)/120;
    third = (float)(z1+z2+2*z3+z4)/120;
    fourth = (float)(z1+z2+z3+2*z4)/120;

    m.at(0).at(0) = first;   m.at(0).at(1) = first;   m.at(0).at(2) = first;
	m.at(1).at(0) = second;   m.at(1).at(1) = second;   m.at(1).at(2) = second; 
    m.at(2).at(0) = third;   m.at(2).at(1) = third;   m.at(2).at(2) = third;
	m.at(3).at(0) = fourth;   m.at(3).at(1) = fourth;   m.at(3).at(2) = fourth; 
    
}

/*CALCULOS PARA PI*/
void calculateGammaPrima(Matrix &m){
    zeroes(m, 12, 3);
    m.at(0).at(0) = (float)1/24;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = (float)1/24;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = (float)1/24;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = (float)1/24;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = (float)1/24;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = (float)1/24;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = (float)1/24;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = (float)1/24;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = (float)1/24;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = (float)1/24; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = (float)1/24;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = (float)1/24; 
}

/*Calculos para Xi*/
float calculateXiFactor(mesh msh, int eIndex){
    element e = msh.getElement(eIndex);

    node node1=msh.getNode(e.getNode1()-1);
    node node2=msh.getNode(e.getNode2()-1);
    node node3=msh.getNode(e.getNode3()-1);
    node node4=msh.getNode(e.getNode4()-1);
    
    float z1, z2, z3, z4;
    
    z1= node1.getZ();
    z2= node2.getZ();
    z3= node3.getZ();
    z4= node4.getZ();
    return (float)pow(z1, 2)+z1*(z2+z3+z4)+pow(z2,2)+z2*(z3+z4)+pow(z3, 2)+z3*z4+pow(z4, 2);
}

float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}
//Casi seguro que check
Matrix createLocalM(int e,mesh &m){
    Matrix matrixChi,matrixZeta,matrixTheta,matrixPi, matrixPI, matrixTao, matrixXi;
    float u_bar,nu,rho,Ve,J,Determinant;
    
    /* [ CHI+zeta  theta + pi ]
       [  PI   tao-Xi ]
    */

    //Matrix Chi
    Matrix Gamma, Alpha, Beta;

    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);

    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    float real_Chi = (float) (J)/(Determinant);
    
    calculateGamma(Gamma, m, e);
    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);
    productRealMatrix(real_Chi, productMatrixMatrix(Gamma,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixChi);

    //Matrix Zeta
    Matrix Alpha_t,Beta_t;
    float real_zeta = (float) (calculateE(m, e)*J)/(Determinant*Determinant);

    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);

    productRealMatrix(real_zeta,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,12),3,3,12),12,3,12),matrixZeta);

    //Matrix theta
    Matrix  A, Q;
    
    float real_theta = (float) (J/(Determinant));
    calculateA(A, m, e);
    calculateQ(Q, m, e);
    productRealMatrix(real_theta,productMatrixMatrix(Q,productMatrixMatrix(Alpha,A,3,3,4),12,3,4),matrixTheta);
    //Matrix pi
    Matrix H;
    float real_Pi = (float)(J/(Determinant));
    calculateH(H, m, e);
    productRealMatrix(real_Pi,productMatrixMatrix(H,productMatrixMatrix(Alpha_t,A,3,3,4),12,3,4),matrixPi);
    //Matrix PI
    Matrix A_t, GammaPrima, GammaPrima_t;
    float real_PI = (float)(J/Determinant);
    calculateGammaPrima(GammaPrima);
    transpose(A, A_t);
    transpose(GammaPrima, GammaPrima_t);
    productRealMatrix(real_PI, productMatrixMatrix(A_t,productMatrixMatrix(Alpha_t, GammaPrima_t, 3, 3, 12),4, 3, 12), matrixPI);
    //Matrix Tao
    Matrix O;
    calculateO(O, m, e);
    float realTao = (float)(J/Determinant);
    productRealMatrix(realTao, productMatrixMatrix(O, productMatrixMatrix(Alpha, A, 3, 3, 4), 4, 3, 4), matrixTao);
    //Matrix Xi
    cout<<calculateXiFactor(m, e)*J<<endl<<pow(Determinant, 1)<<endl;
    float real_xi = (float)(calculateXiFactor(m, e)*J/30*(Determinant*Determinant));
    cout<<real_xi<<endl;
    productRealMatrix(real_xi, productMatrixMatrix(A_t, productMatrixMatrix(Alpha_t, productMatrixMatrix(Alpha, A, 3, 3, 4), 3, 3, 4), 4, 3, 4), matrixXi);
    //Matrix M
    Matrix M;
    zeroes(M,16);
    ubicarSubMatriz(M,0,11,0,11, sumMatrix(matrixChi,matrixZeta,12,12));
    ubicarSubMatriz(M,0, 11, 12, 15, sumMatrix(matrixTheta, matrixPi, 12, 4));
    ubicarSubMatriz(M, 12, 15, 0, 11, matrixPI);
    ubicarSubMatriz(M,12,15,12,15,resMatrix(matrixTheta, matrixXi, 4, 4));
    return M;
}

void calculateGravity(Vector &f, mesh &m){
    zeroes(f,3);
    f.at(0) += m.getParameter(EXTERNAL_FORCE_X);
    f.at(1) += m.getParameter(EXTERNAL_FORCE_Y);
    f.at(2) += m.getParameter(EXTERNAL_FORCE_Z);

}

void calculateL(Vector &f, mesh m, int eIndex){
    zeroes(f, 4);
    element e = m.getElement(eIndex);
    node node1=m.getNode(e.getNode1()-1);
    node node2=m.getNode(e.getNode2()-1);
    node node3=m.getNode(e.getNode3()-1);
    node node4=m.getNode(e.getNode4()-1);
    
    float z1, z2, z3, z4;
    
    z1= node1.getY();
    z2= node2.getY();
    z3= node3.getY();
    z4= node4.getY();
    float first, second, third, fourth;

    first = (float)7*(2*z1+z2+z3+z4)/120;
    second = (float)7*(z1+2*z2+z3+z4)/120;
    third =(float) 7*(z1+z2+2*z3+z4)/120;
    fourth = (float)7*(z1+z2+z3+2*z4)/120;
    f.at(0)=first;f.at(1)=second;f.at(2)=third; f.at(3)=fourth;
}

Vector createLocalb(int e,mesh &m){
    float J;
    Vector b,b_aux,gravity, zvector;
    Matrix GammaPrima;

    calculateGravity(gravity, m);
    calculateL(zvector, m, e);
    calculateGammaPrima(GammaPrima);
    

    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    zeroes(b_aux,16);
    productMatrixVector(GammaPrima,gravity,b_aux);
    b_aux.at(12) = zvector.at(0);
    b_aux.at(13) = zvector.at(1);
    b_aux.at(14) = zvector.at(2);
    b_aux.at(15) = zvector.at(3);
    productRealVector(J,b_aux,b);
    return b;
}
