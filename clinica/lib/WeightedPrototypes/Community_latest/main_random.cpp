using namespace std;

char *outfile  = NULL;

int
main(int argc, char **argv) {
  srand(time(NULL)+getpid());

  int n = atoi(argv[1]);
  int degree = atoi(argv[2]);

  for (unsigned int i=0 ; i<n ; i++) {
    for (unsigned int j=0 ; j<n ; j++) {
      int r  = rand()%n;
      if (r<degree)
        cout << i << " " << j << endl;
    }
  }
}
