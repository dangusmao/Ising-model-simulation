#include<stdio.h>
#include<math.h>
#include<unistd.h>
#include<stdlib.h>
#include<time.h>

							//Modelo de Ising 3-D - Simulação de Monte Carlo - DEBUG

									//Aluno: Daniel de Gusmão Moraes

//Assumindo que a termalização ocorre após 100.000 passos, tiramos as médias apenas nos 1000 últimos passos.

			#define P 110000

//Tamanho do vetor usado pra tirar as médias

			#define TAM 10000 // 10 mil

//Dimensão da matriz 3-D

			#define M 30
													//||\\

FILE* saidaE[20]; //Saída das energias de cada Beta
FILE* saidaM[20]; //Saída das magnetizações de cada Beta
FILE* sE;		  //Saída dos valores médios de energia
FILE* sM;		  //Saída dos valores médios de magnetização


float avg(float A[TAM]){

	int i; float a=0;

	for (i = 0; i < TAM; ++i)
	{
		a+=A[i];
	}

	return a/(float)TAM;
}

float sd(float A[TAM],float avg){

	int i; float a=0;

	for (i = 0; i < TAM; ++i)
	{
		a+=pow(A[i]-avg,2);
	}

	a=a/(float)TAM;	a=sqrt(a);

	return a;

}

int sortspin(){
	
	int s;

	s = ( rand() % (M) ) + 0;

	return s;
}

float energia(int lat[M][M][M]){
	
	int i,j,k;	float h=0;

	for (i = 0; i < M; ++i)
	{
		for (j = 0; j < M; ++j)
		{
			for (k = 0; k < M; ++k)
			{
				h+=lat[i][j][k]*lat[(i+1)%M][j][k];
				h+=lat[i][j][k]*lat[i][(j+1)%M][k];
				h+=lat[i][j][k]*lat[i][j][(k+1)%M];
			}
		}
	}

	return -h;

}

float mag(int lat[M][M][M]){

	int i,j,k;	float m=0;

	for (i = 0; i < M; ++i)
	{
		for (j = 0; j < M; ++j)
		{
			for (k = 0; k < M; ++k)
			{
				m+=(float)lat[i][j][k];
			}			
		}
	}

	m=m/((float)M*(float)M*(float)M);

	return m;
}

float dE(int si,int sj,int sk,int lat[M][M][M]){

	float dE;

	dE = lat[si][sj][sk]*lat[(si+1)%M][sj][sk] + lat[si][sj][sk]*lat[si][(sj+1)%M][sk] + lat[si][sj][sk]*lat[si][sj][(sk+1)%M];

	if (si==0){	dE+=lat[M-1][sj][sk]*lat[si][sj][sk];} else dE+= lat[(si-1)%M][sj][sk]  *lat[si][sj][sk];
	if (sj==0){	dE+=lat[si][M-1][sk]*lat[si][sj][sk];} else dE+= lat[si][sj-1][sk]      *lat[si][sj][sk];
	if (sk==0){	dE+=lat[si][sj][M-1]*lat[si][sj][sk];} else dE+= lat[si][sj][(si-1)%M]  *lat[si][sj][sk];

//VERIFICAR NO DEBUG: TALVEZ SEJA -2dE mesmo

	dE=-2*dE; return dE;

}



void metropolis(int f,float B, int lat[M][M][M], float Eint[], float Mag[]){

	int i,si,sj,sk;		float Eo,E,deltaE; 		double p,bp; 

	Eo = energia(lat); 

	for (i = 0; i < P; ++i)
	{		
		//Impressão de energias e mag. individuais

		fprintf(saidaE[f],"%d\t%f\n",i,Eo);
		fprintf(saidaM[f],"%d\t%f\n",i,mag(lat));

		//Verifica se o equiíbrio já foi atingido

		if (i>=P-TAM-1)
		{
			Eint[i-(P-TAM-1)] = Eo ; Mag[i-(P-TAM-1)] = mag(lat);
		}
		
		//Passo de Monte Carlo		

		si = sortspin();	sj = sortspin();	sk = sortspin();	lat[si][sj][sk]=-1*lat[si][sj][sk];			

		deltaE = dE(si,sj,sk,lat);		bp = exp(-1*B*deltaE);		
		

		if (deltaE<0)
		{
			Eo=Eo+deltaE;		
		}

		if (deltaE>0)
		{
			p=(float)rand() / (float)RAND_MAX;
			
			if (p<=bp)
			{
				Eo=Eo+deltaE;							
			}

			if (p>bp)
			{
				lat[si][sj][sk]=-lat[si][sj][sk];
								
			}
		}
	}
}


int main(){	


	sE = fopen("mediaEnergiaBeta","w");
	sM = fopen("mediaMagBeta","w");  

	srand(time(NULL)); 

	int i,j,k,c,lat[M][M][M];		float E[20000],Mag[20000]; 	float B;	 float avgE,sdE,avgM,sdM;

	char cmdE[100];
	char cmdM[100];	//vetores criado para facilitar a mudança dos nomes dos arquivos de saída	

	//Simulações para cada Beta

	B=0.00;
	
	for(c=0;c<21;c++)
	{	
		//Inicialização dos spins

		for (i = 0; i < M; ++i)
		{
			for (j = 0; j < M; ++j)
			{
				for (k = 0; k < M; ++k)
				{
					lat[i][j][k]=1;
				}				
			}
		}
				
		saidaE[c] = fopen("Energia.dat","w");
		saidaM[c] = fopen("Magnet.dat","w");	

		metropolis(c,B,lat,E,Mag);

		//Médias + Desvio padrão

		avgE=avg(E); avgM=avg(Mag) ; sdE=sd(E,avgE) ; sdM=sd(Mag,avgM);

		fprintf(sE,"%f\t%f\t%f\n",B,avgE,sdE);
		fprintf(sM,"%f\t%f\t%f\n",B,avgM,sdM);			

		sprintf(cmdE,"mv Energia.dat Energia%d.dat",c);		sprintf(cmdM,"mv Magnet.dat Magnet%d.dat",c);
		
		fclose(saidaE[c]);	fclose(saidaM[c]);		
				
		system(cmdE);	system(cmdM);

		printf("%d/20\n",c+1);

		B+=0.05;
		
	}

	fclose(sE);fclose(sM);

	return 0;
}

