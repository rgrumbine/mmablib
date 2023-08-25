#include "ncepgrids.h"

// Print out grid information a la 
//     https://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID172

void show(ijpt &fij, latpt &ll) ;

int main(void) {
	GRIDTYPE<float> x;
	ijpt fij;
	latpt ll;

	printf("ni = %3d\n",x.xpoints() );
	printf("nj = %3d\n",x.ypoints() );

	fij.i = 0;
	fij.j = 0;
        ll = x.locate(fij);
	show(fij, ll);

	fij.j = x.ypoints()- 1;
        ll = x.locate(fij);
	show(fij, ll);

	fij.i = x.xpoints() - 1;
	ll = x.locate(fij);
	show(fij, ll);

	fij.j = 0;
	ll = x.locate(fij);
	show(fij, ll);

	ll.lat = -90.0;
	ll.lon = 0.0;
	fijpt ffij;
  	ffij = x.locate(ll);
  	printf("south pole is %f %f\n",ffij.i, ffij.j);

	int delta = 40;
	for (fij.i = x.xpoints()/2-delta; fij.i <= x.xpoints()/2+delta ; fij.i++) {
	for (fij.j = x.ypoints()/2-delta; fij.j <= x.ypoints()/2+delta ; fij.j++) {
	  ll = x.locate(fij); 
	  if (ll.lat < -89.75) {
	    printf("%3d %3d pt is: %f %f\n",(int)fij.i+1, (int)fij.j+1, ll.lon, ll.lat);
	  }
	}
	}
	
	return 0;
}
void show(ijpt &fij, latpt &ll) {
  char pole = 'N';
  char hemi = 'E';

  if (ll.lat < 0.0) {
    pole = (char) 'S';
  }
  if (ll.lon > 180.0) {
    ll.lon -= 360.;
  }
  else if (ll.lon < -180.0) {
    ll.lon += 360.;
  }

  if (ll.lon < 0.) {
    hemi = (char) 'W';
  }

  printf("i,j pt: (%3d,%3d) lon,lat: %9.5f %1c %9.5f %1c\n",(int)fij.i+1, (int)fij.j+1, 
	  fabs(ll.lon),hemi, fabs(ll.lat),pole);

  return;
}
