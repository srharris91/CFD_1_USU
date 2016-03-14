// use point SOR to solve the finite-volume method "Del^2 phi = 0"
// Ax=b
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class node{
    public:
        double x,y,n,phi,phiN;  // x,y coordinates,node number, and phi value, and next iteration phi value
        double N,E,W,S; // neighbors in North, East, West, and South directions for 2d
        node(){             // constructor of 0's
            x       =   0.  ;
            y       =   0.  ;
            n       =   0.  ;
            phi     =   0.  ;
            phiN    =   0.  ;
        }


};
void show(vector< vector<node> > a){
    cout<<" x,y,n,phi,phiN for each position = "<<endl;
    for (auto nodes_it = a.begin(); nodes_it!=a.end(); ++nodes_it){
        for (auto  nodes_it_j = nodes_it->begin(); nodes_it_j!=nodes_it->end(); ++nodes_it_j){
            cout
                <<nodes_it_j->x<<","
                <<nodes_it_j->y<<","
                <<nodes_it_j->n<<","
                <<nodes_it_j->phi<<","
                <<nodes_it_j->phiN<<","
                <<"  ";
        }
        cout<<"end line"<<endl;
    }
}

int main(){
    vector< vector<node> > nodes(20,vector<node> (20));       // 20 by 20 matrix containing node data
    vector< vector<double> > A,D,L,U;   // A matrix and it's subcomponents D=diagonal,L=lower,U=upper (A=D + L + U)
    vector<double> x,b;
    vector<double> aw,ae,as,an,ap;         // coefficients on west, east, south and north neighbors of center point
    vector< vector<double> >::iterator i;
    vector<double>::iterator j;
    double temp;
    double w;       // relaxation factor for SOR method
    double dx=1./20.,dy=1./20.;     // spacing

    int MAX_R=20,MAX_C=20;    // 20 by 20 finite volume 2d problem

    // set node values
    int xi=0, yj=0, nn=0;
    for (auto nodes_it = nodes.begin(); nodes_it!=nodes.end(); ++nodes_it,++yj){
        xi=0;
        for (auto  nodes_it_j = nodes_it->begin(); nodes_it_j!=nodes_it->end(); ++nodes_it_j,++xi){
            nodes_it_j->n   = nn++;
            nodes_it_j->phi = 0.;
            nodes_it_j->x   = xi*dx+dx/2.;
            nodes_it_j->y   = yj*dy+dy/2.;
            if (xi*dx==1.) nodes_it_j->phi = yj;
            if (yj*dx==1.) nodes_it_j->phi = xi;
        }
    }
    // show nodes
    show(nodes);
    // set aw, ae, as, and an, and ap values 
    ap.push_back(2.+0.);    // top left BC
    for (int xi=0;xi<MAX_R; xi++){
        for (int yj=0;yj<MAX_C; yj++){
            aw.push_back(1.);
            ae.push_back(1.);
            as.push_back(1.);
            an.push_back(1.);
            ap.push_back(aw[xi]+ae[xi]+as[xi]+an[xi]);
        }
    }
    // set BC values on RHS and diagonal(ap)



    //initialize values of A (D,L, and U subcomponents) and b
    /*        for (int xi=0; xi<MAX_R; xi++){
              A.push_back(vector<double> ());
              for (int yj=0; yj<MAX_C; yj++){
              A[xi].push_back(double(xi*yj));
    //            cout<<A[xi][yj]<<" ";
    }
    //        cout<<endl;
    }
    //    cout<<endl<<endl;
    // split into subcomponents
    xi=0, yj=0;
    for (i=A.begin(); i!=A.end(); i++,xi++){
    D.push_back(vector<double> ());
    U.push_back(vector<double> ());
    L.push_back(vector<double> ());
    yj = 0;
    for (j=i->begin(); j!=i->end(); j++,yj++){
    temp = *j;
    if (xi == yj){              // diagonal
    D[xi].push_back(temp);
    U[xi].push_back(0.);
    L[xi].push_back(0.);
    }
    else if (xi < yj) {         // upper
    D[xi].push_back(0.);
    U[xi].push_back(*j);
    L[xi].push_back(0.);
    }
    else if (xi > yj){          // lower
    D[xi].push_back(0.);
    U[xi].push_back(0.);
    L[xi].push_back(*j);
    }
    else 
    cout<<"something is wrong"<<endl;

    //            cout<<*j<<" ";
    }
    //        cout<<endl;
    }
    // output to check
    //    cout<<"D = "<<endl;
    //    for (xi=0;xi<MAX_R;++xi){
    //        for (yj=0;yj<MAX_C;++yj){
    //            cout<<D[xi][yj]<<" ";
    //        }
    //        cout<<endl;
    //    }
    //    cout<<"U = "<<endl;
    //    for (xi=0;xi<MAX_R;++xi){
    //        for (yj=0;yj<MAX_C;++yj){
    //            cout<<U[xi][yj]<<" ";
    //        }
    //        cout<<endl;
    //    }
    //    cout<<"L = "<<endl;
    //    for (xi=0;xi<MAX_R;++xi){
    //        for (yj=0;yj<MAX_C;++yj){
    //            cout<<L[xi][yj]<<" ";
    //        }
    //        cout<<endl;
    //    }



    // solve x values using SOR method


*/
    return 0;
}
