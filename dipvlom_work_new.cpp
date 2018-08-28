// Построение каталога гравитационных моделей для различных астероидов

#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>

using namespace std;
typedef vector<double> vecd;
typedef vector<vecd> state_type;

#define EXP 2.7182818284590
#define PI 3.141592653589

class m_point {
private:
    double r_,theta_,lambda_;
    bool sph_r, sph_theta, sph_lambda;
public:
    vecd coord;
    double kappa;
    m_point(double x, double y, double z): coord{x,y,z},kappa(1),sph_r(false), sph_theta(false), sph_lambda(false) {}
    double calc_pot(const vecd& prob){
        double dist=0;
        for (int j=0;j<prob.size();j++)
            dist += (prob[j] - coord[j])*(prob[j] - coord[j]);
        dist = sqrt(dist);
        return kappa / dist;
    }
    double r(){
        if (!sph_r){
            r_ = sqrt(coord[0]*coord[0]+coord[1]*coord[1]+coord[2]*coord[2]);
            sph_r = true;
        }
        return r_;
    }
    double theta(){
        if (!sph_theta){
            theta_ = atan2(sqrt(coord[0]*coord[0]+coord[1]*coord[1]),coord[2]);
            sph_theta = true;
        }
        return theta_;
    }
    double lambda(){
        if (!sph_lambda){
            lambda_ = atan2(coord[1],coord[0]);
            sph_lambda = true;
        }
        return lambda_;
    }
};

struct vertex {
    double x, y, z;
    double dist() {
        return max({abs(x), abs(y), abs(z)});
    }
    double rms() {
        return sqrt(x*x+y*y+z*z);
    }
};

class poligon {
private:
    vertex normal, ver;
public:
    poligon(vertex v1, vertex v2, vertex v3): ver(v1){
        vertex a, b;
        a.x = v2.x - v1.x;
        a.y = v2.y - v1.y;
        a.z = v2.z - v1.z;
        b.x = v3.x - v2.x;
        b.y = v3.y - v2.y;
        b.z = v3.z - v2.z;
        normal.x = a.y*b.z - a.z*b.y;
        normal.y = a.z*b.x - a.x*b.z;
        normal.z = a.x*b.y - a.y*b.x;
    }
    bool right_side(double x, double y, double z){
        return 0 >= normal.x * (x - ver.x) + normal.y * (y - ver.y) + normal.z * (z - ver.z);
    }
};

class body{
private:
    vector<poligon> poligons;
public:
    double maximum;
    double scale_factor;
    void open(const string& str){
        ifstream aster_data(str);
        vector<vertex> vertices{};
        char l;
        aster_data >> l;
        while (l == 'v') {
            vertex v{};
            aster_data >> v.x >> v.y >> v.z;
            vertices.emplace_back(v);
            aster_data >> l;
        }
        cout << "number of vertices " << vertices.size() << endl;
        while (l == 'f') {
            int i, j, k;
            aster_data >> i >> j >> k;
            poligons.emplace_back(poligon(vertices[i-1], vertices[j-1], vertices[k-1]));
            aster_data >> l;
        }
        cout << "number of poligons " << poligons.size() << endl;
        aster_data.close();
        for (auto& v:vertices) {
            if (v.dist() > maximum)
                maximum = v.dist();
            if (v.rms() > scale_factor)
                scale_factor = v.rms();
        }
    }

    bool is_inside(double x, double y, double z){
        for (auto& p:poligons)
            if (!p.right_side(x, y, z))
                return false;
        return true;
    }

};



//class ellipsoid {
//    double a_, b_, c_;
//public:
//    ellipsoid(double a,double b,double c): a_(a), b_(b), c_(c){}
//    bool is_inside(double x, double y, double z) {
//        return 0 >= x*x/(a_*a_) + y*y/(b_*b_) + z*z/(c_*c_) - 1;
////        return 0 >= (x*x+y*y+z*z)*(x*x+y*y+z*z) - 2*c_*c_*(x*x-y*y-z*z) - 2*pow(a_,4) + pow(c_,4);
//    }
//};

double kron_delta(double n, double k){
    if (n == k)
        return 1;
    else
        return 0;
}

double sgn(double n){
    if (n < 0)
        return -1;
    else
        if (n > 0)
            return 1;
        else
            return 0;
}

state_type mult_mat(const state_type& A, const state_type& B){
    assert(A[0].size() != B.size());
    state_type C(A.size());
    for (auto& str:C)
        str.resize(B[0].size());
    for (int i=0; i<A.size();i++)
        for (int j=0;j<B[0].size();j++)
            for (int k=0;k<A[0].size();k++)
                C[i][j] += A[i][k]*B[k][j];
    return C;
}

state_type transpose_mat(const state_type& A){
    state_type C(A[0].size());
    for (auto& str:C)
        str.resize(A.size());
    for (int i=0; i<A[0].size();i++)
        for (int j=0;j<A.size();j++)
            C[i][j] = A[j][i];
    return C;
}

state_type matrix_rot(int elem_i, int elem_j, int dim, double theta){
    state_type A(dim);
    for (auto& str:A)
        str.resize(dim);
    for (int i=0;i<dim;i++)
        for (int j=0;j<dim;j++)
            if (i == j)
                A[i][j] = 1;
            else A[i][j] = 0;
    A[elem_i][elem_i] = A[elem_j][elem_j] = cos(theta);
    A[elem_j][elem_i] = -A[elem_i][elem_j] = -sin(theta);
};

class tensor{
    state_type matrix;
    void rot_m(state_type& m){
        double eps = numeric_limits<double>::epsilon();
        double t_m = t_diag = 1;
        double tau,t,c,s;
        while (t_m/t_diag > eps){
            double max_el=0, i_max, j_max;
            t_m = t_diag = 0;
            for (int i=0;i<m.size();i++)
                for (int j=0;j<m.size();j++)
                    if (i != j){
                        if (abs(m[i][j]) > max_el){
                            max_el = abs(m[i][j]);
                            i_max = i;
                            j_max = j;
                        }
                        t_m += m[i][j];
                    }
                    else
                        t_diag += m[i][j];
            if (m[i_max][i_max] == m[j_max][j_max])
                theta = PI/4;
            else {
                tau = (a(j,j)-a(k,k)/2/a(j,k)
                t = sign(tau)/(|tau|+sqrt(1+tau^2))
                c = 1/sqrt(1+t^2)
                s = ct
            }

        }
    }
public:
    tensor(const vector<m_point>& mps): matrix(3) {
        for (auto& ten:matrix)
            ten.resize(3);
        for (auto mp:mps){
            matrix[0][0] += mp.coord[1]*mp.coord[1] + mp.coord[2]*mp.coord[2];
            matrix[1][1] += mp.coord[0]*mp.coord[0] + mp.coord[2]*mp.coord[2];
            matrix[2][2] += mp.coord[0]*mp.coord[0] + mp.coord[1]*mp.coord[1];
            matrix[0][1] = matrix[1][0] -= mp.coord[0] * mp.coord[1];
            matrix[0][2] = matrix[2][0] -= mp.coord[0] * mp.coord[2];
            matrix[1][2] = matrix[2][1] -= mp.coord[1] * mp.coord[2];
        }
        for (auto& ten:matrix)
            for (auto& elem:ten)
                elem *= mps.size();
    }

    void transform(vector<m_point>& mps){
        state_type m(3);
        for(auto& ten:m)
            ten.resize(3);
        rot_m(m);
        for (auto& mp:mps){
            m_point copy_mp(mp);
            for (int i=0;i<3;i++){
                mp.coord[i] = 0;
                for (int j=0;j<3;j++)
                    mp.coord[i] += m[i][j] * copy_mp[j];
            }
        }
    }
};

double Leg_pol(double n, double k, double phi) {
    double cos_phi=cos(phi), sin_phi=sin(phi);
    double pol_prev=0, pol_next, pol=1;

    for (int j=0;j<k;j++) {
        pol *= -(2*j+1)*abs(sin_phi)*sqrt((1+kron_delta(0,j))*static_cast<double>(2*j+3)/(2*j+1)/(2*j+1)/(2*j+2));
    }
    for (int i=k;i<n;i++) {
        pol_next = ((2*i+1)*cos_phi*pol*sqrt(static_cast<double>(2*i+3)/(2*i+1)*(i-k+1)/(i+k+1)) - (pol_prev == 0 ? 0:(i+k)*pol_prev*sqrt(static_cast<double>(2*i+3)/(2*i-1)*(i-k+1)/(i+k+1)*(i-k)/(i+k))))/(i-k+1);
        pol_prev = pol;
        pol = pol_next;
    }
    return pol;
}

double inv_alpha(double n, double k){
    double mult=1;
    for (int i=n-k+1;i<=n+k;i++)
        mult*=i;
    return mult/(2-kron_delta(0,k));
}

//double fact(int N)
//{
//    if(N < 0)
//        return 0;
//    if (N == 0)
//        return 1;
//    else
//        return N * fact(N - 1);
//}

//double det_max (double a, double b, double c) {
//    double d;
//    d = a;
//    if (a < b) {
//        d = b;
//        if (b < c) {
//            d = c;
//        }
//    } else if (a < c) {
//        d = c;
//    }
//    return d;
//}

//double funy(double x, double a, double b){
//    double d;
//    d = sqrt(b*b - b*b/(a*a)*x*x);
//    return d;
//}
//
//double funz(double x, double y, double a, double b, double c){
//    double e;
//    e = sqrt(c*c - (c*c/(a*a)*x*x + c*c/(b*b)*y*y));
//    return e;
//}

int main(){
    double x, y, z, h, sum_pot=0, r;
    int n, num_coef;

//    cout << Leg_pol(10,0,0.6) << endl;
//    exit(0);

//    vecd gradU(3);
    std::ofstream out, res, graph, file_coef;
    std::ifstream file_in, aster_data;

    file_in.open("file_in.dat",std::ios_base::in);
//  ввод параметров эллипса и шага
    file_in >> n >> num_coef;
    file_in.close();

    out.open("out.dat",std::ios_base::out);

//    vector<ellipsoid> ell_body{};

//    file_coef.open("file_coef.dat",std::ios_base::in);
//    for (int i=0;i<k;i++) {
//        file_coef >> a >> b >> c;
//        ell_body.emplace_back(ellipsoid(a, b, c));
//        if (max({a,b,c}) > m) {
//            m = max({a,b,c});
//        }
//    }
//    file_coef.close();

//    aster_data.open("aster_data.dat",std::ios_base::in);
//    vector<vertex> vertices{};
//    char l;
//    aster_data >> l;
//    while (l == 'v') {
//        vertex v{};
//        aster_data >> v.x >> v.y >> v.z;
//        vertices.emplace_back(v);
//        aster_data >> l;
//    }
//    cout << "number of vertices " << vertices.size() << endl;
//    vector<poligon> poligons{};
//    while (l == 'f') {
//        int i, j, k;
//        aster_data >> i >> j >> k;
//        poligons.emplace_back(poligon(vertices[i], vertices[j], vertices[k]));
//        aster_data >> l;
//    }
//    cout << "number of poligons " << poligons.size() << endl;
//    aster_data.close();

    body aster{};
    aster.open("aster_data.dat");

    double m = aster.maximum;
    h = 2.0 * m / (n-1);

    vector<m_point> disc_body{};
    disc_body.reserve(n*n*n);

    for (x=-m; x<=m; x+=h) {
        for (y=-m; y<=m; y+=h) {
            for (z=-m; z<=m; z+=h) {
                if (aster.is_inside(x, y, z)) {
//                    out << x << " " << y << " " << z << endl;
                    disc_body.emplace_back(m_point(x,y,z));

                }
            }
        }
    }

    vecd center_m(3);
    for (auto& p:disc_body){
        for (int i=0;i<3;i++){
            center_m[i] += p.coord[i]/disc_body.size();
        }
    }
//    cout << center_m[0] << " " << center_m[1] << " " << center_m[2] << endl;
//    out << endl << "new" << endl;
    for (auto& p:disc_body){
        for (int i=0;i<3;i++){
            p.coord[i] -= center_m[i];
            out << p.coord[i] << " ";
        }
        out << endl;
    }



    state_type A(num_coef), B(num_coef);
    file_coef.open("file_coef.dat",std::ios_base::out);
    for (int i=0;i<=num_coef;i++) {
        A[i].resize(i+1);
        B[i].resize(i+1);
        for (int j=0;j<=i;j++){
            A[i][j] = 0;
            B[i][j] = 0;
            for (auto& p:disc_body){
                A[i][j] += pow(p.r(),i)*Leg_pol(i,j,p.theta())*cos(j*p.lambda());
                B[i][j] += pow(p.r(),i)*Leg_pol(i,j,p.theta())*sin(j*p.lambda());
            }
            A[i][j] /= inv_alpha(i,j)*pow(aster.scale_factor,i)*disc_body.size();
            B[i][j] /= inv_alpha(i,j)*pow(aster.scale_factor,i)*disc_body.size();
            file_coef << i << " " << j << " " << A[i][j] << " " << B[i][j] << endl;
            cout << A[i][j] << " " << B[i][j] << " ";
        }
        cout << endl;
    }
    file_coef.close();

    res.open("res.dat", std::ios_base::out);
    graph.open("graph.dat",std::ios_base::out);

    r = 300;
    vecd prob_p(3);
    prob_p[0] = 0; prob_p[1] = 0; prob_p[2] = -r;
    for (int i=0;i<disc_body.size();i++) {
        sum_pot += disc_body[i].calc_pot(prob_p);
    }
    res << "theta= " << 180 << " " << "lambda= " << 0 << " sum_pot= " << sum_pot << endl;
    for (int theta=179;theta>0;theta-=1) {
        for (int lambda=0; lambda<=360; lambda+=1) {
            sum_pot=0;
            prob_p[0]=r*sin((double) theta*PI/180.0)*cos((double) lambda*PI/180.0);
            prob_p[1]=r*sin((double) theta*PI/180.0)*sin((double) lambda*PI/180.0);
            prob_p[2]=r*cos((double) theta*PI/180.0);
            for (int i=0;i<disc_body.size();i++) {
                sum_pot += disc_body[i].calc_pot(prob_p);
            }
            res << "theta= " << theta << " lambda= " << lambda << " sum_pot= " << sum_pot << endl;
            graph << theta << " " << lambda << " " << sum_pot << endl;
        }
    }
    sum_pot=0;
    prob_p[0] = 0; prob_p[1] = 0; prob_p[2] = r;
    for (int i=0;i<disc_body.size();i++) {
        sum_pot += disc_body[i].calc_pot(prob_p);
    }
    res << "theta= " << 0 << " lambda= " << 0 << " sum_pot= " << sum_pot << endl;

//    while (l.is_inside(x,y,z))
//        while (l.is_inside(x,y,z)){
//            while (l.is_inside(x,y,z)) {
//                out << x << " " << y << " " << z << endl;
//                x+=h;
//            }
//            y+=h;
//        }
//        z+=h;
//    }

//    for (x=-m; x<=m; x+=h) {
//        for (y=-m; y<=m; y+=h) {
//            for (z=-m; z<=m; z+=h) {
//                if (any_of(ell_body.begin(),ell_body.end(),[=](ellipsoid ell){return ell.is_inside(x,y,z);})) {
//                out << x << " " << y << " " << z << endl;
//                disc_body.emplace_back(m_point(x,y,z));
//                }
//            }
//        }
//    }

    out.close();
    res.close();
    graph.close();

    return 0;
}


