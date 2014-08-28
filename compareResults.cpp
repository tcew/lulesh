
#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>

int main(int argc, char *argv[]){

  if(argc != 4){
    printf("Usage is ./compareResults n <filename1> <filename2> \n");
    exit(-1);
  }

  const int n = atoi(argv[1]);
  const int N = (n+1)*(n+1)*(n+1);

  std::string f1 = argv[2];
  std::string f2 = argv[3];

  std::vector<double> v1(3*N, 0.);
  std::vector<double> v2(3*N, 0.);

  std::ifstream file;

  file.open(f1.c_str(), std::ifstream::in);

  for(int i=0; i<N; i++){
    int dummy;

    file >> dummy >> v1[3*i] >> v1[3*i+1] >> v1[3*i+2];

  }

  file.close();

  file.open(f2.c_str(), std::ifstream::in);

  for(int i=0; i<N; i++){
    int dummy;

    file >> dummy >> v2[3*i] >> v2[3*i+1] >> v2[3*i+2];

  }

  file.close();

  double lInfError = 0.;
  double l2Error = 0.;

  for(int i=0; i<3*N; i++){

    double err = fabs(v1[i]-v2[i]);

    l2Error += err*err;

    if(lInfError < err)
      lInfError = err;
  }

  printf("Error stats:\n");
  printf("\t l2Error   = %lg \n", sqrt(l2Error));
  printf("\t lInfError = %lg \n", lInfError);

  return 0;

}
