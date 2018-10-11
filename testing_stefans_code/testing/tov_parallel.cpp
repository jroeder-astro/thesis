#include <iostream>                 //Ein-/Ausgabe (Include-Dateien)
#include <math.h>                   //Mathematisches
#include <omp.h>                    //OpenMP
#include <stdio.h>                  //Fuer die Ausgabedatei

//Definition der Zustandsgleichung
double eos(double p)
  {
    double e;
    e=pow(p/10,3.0/5);
    return e;
  }  

main(void)  //Hauptprogramm
  {
    //Variablendeklarationen
    int i,anz=40;
    double M,p,e,r,nu,dM,dp,de,dr,dnu,dec;
    double Er[anz],EM[anz],Enu[anz],Eec[anz];
    double eos(double);
    
    //Ausgabedatei
    FILE *ausgabe;
    ausgabe = fopen("output/tov.txt", "w+");           
    
    //Variableninitialisierung
    dr=0.000001;
    dec=0.0001;
    
    //for Schleife zur Berechnung mehrerer Sterne

#pragma omp parallel for private(i,M,p,e,r,nu,dM,dp,de,dnu)

    for (i=0;i<40;i++)
    {
      M=0;
      r=pow(10,-14);
      p=10*pow(0.0005+i*dec,5.0/3);
      nu=0;
      Eec[i]=eos(p);
      
      //do-while Schleife (Numerische Lösung der TOV-Gleichung) 
      do
        {
          e=eos(p);                                      
          dM=4*M_PI*e*r*r*dr;                            
          dp=-(p+e)*(M+4*M_PI*r*r*r*p)/(r*(r-2*M))*dr;   
          dnu=(M + 4*M_PI*r*r*r*p)/(r*(r-2*M))*dr;       
          r=r+dr;                                        
          M=M+dM;                                        
          p=p+dp;                                        
          nu=nu+dnu;                                     
        }
      while(p>0);
      
      Er[i]=r;
      EM[i]=M;
      Enu[i]=log(1-2*M/r)/2-nu; 
    }
    
    //Geordnete Ausgabe der Masse, des Radius und der zentralen g00-Metrikkomponente in die Ausgabedatei
    fprintf(ausgabe,"# R[km]    M[Msol]   g00    ec[MeV/fm3] \n"); 
    for (i=0;i<40;i++)
    {
      fprintf(ausgabe, "%f %f %f %f \n",Er[i],EM[i]/1.4766,exp(2*Enu[i]),Eec[i]*pow(10,6)/1.3234); 
    }
    
    fclose(ausgabe);                                   //Ausgabedatei schliessen
    
    return 0;                                          //main beenden (Programmende)
  }
