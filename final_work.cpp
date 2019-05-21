// Building a catalog of gravity models for asteroids with a known form

#include <iostream>
#include <ctime>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <string>

using namespace std;
typedef vector<double> state_vector;
typedef vector<state_vector> state_matrix;

#define EXP 2.7182818284590
#define PI 3.141592653589

class m_point {
private:
    double r_,theta_,lambda_;
    bool sph_r, sph_theta, sph_lambda;
public:
    state_vector coord;
    double kappa;
    m_point(double x, double y, double z): coord{x,y,z},kappa(1),sph_r(false), sph_theta(false), sph_lambda(false) {}
    double calc_pot(const state_vector& prob){
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
        return max({fabs(x), fabs(y), fabs(z)});
    }
    double rms() {
        return sqrt(x*x+y*y+z*z);
    }
};

class polygon {
private:
    vertex normal, ver1, ver2, ver3;
    double d;
public:
    polygon(vertex v1, vertex v2, vertex v3): ver1(v1), ver2(v2), ver3(v3){
        vertex a, b;
        a.x = v2.x - v1.x;
        a.y = v2.y - v1.y;
        a.z = v2.z - v1.z;
        b.x = v3.x - v2.x;
        b.y = v3.y - v2.y;
        b.z = v3.z - v2.z;
        normal.x = a.y*b.z - a.z*b.y;                       // a
        normal.y = a.z*b.x - a.x*b.z;                       // b
        normal.z = a.x*b.y - a.y*b.x;                       // c
        d = normal.x*v1.x + normal.y*v1.y + normal.z*v1.z;  // d
    }

    double fun_z(double i, double j){
            return (d - normal.x * i - normal.y * j) / normal.z;
    }

    bool right_side(double x0, double y0, double z0){
        double eta, dzeta, psi;
        double k, l, rel;
        eta = (ver1.x - x0)*(ver2.y - ver1.y) - (ver2.x - ver1.x)*(ver1.y - y0);
        dzeta = (ver2.x - x0)*(ver3.y - ver2.y) - (ver3.x - ver2.x)*(ver2.y - y0);
        psi = (ver3.x - x0)*(ver1.y - ver3.y) - (ver1.x - ver3.x)*(ver3.y - y0);
        if ((eta*dzeta >= 0) && (dzeta*psi >= 0) && (eta*psi >= 0)){
//            k = min(0.0, fun_z(x0, y0));
//            l = max(0.0, fun_z(x0, y0));
            //cout << k << " " << l << endl;
            rel = z0/fun_z(x0,y0);
            if ((rel >= 0) && (rel < 1))
                return true;
            else
                return false;
        }
        else
            return false;
    }

//---------------- old --------------------------
//    bool right_side(double x, double y, double z){
//        return 0 >= normal.x * (x - ver1.x) + normal.y * (y - ver1.y) + normal.z * (z - ver1.z);
//    }

    double area(){
        return sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
    }

};

class body{
private:
    vector<polygon> polygons;
public:
    double maximum;
    double scale_factor;
    bool existence = true;
    void open(const string& str, const string& vert){
        ifstream aster_data(str);
        ofstream ver_data(vert);
        vector<vertex> vertices{};
        char l;
        aster_data >> l;
        while (l == 'v') {
            vertex v{};
            aster_data >> v.x >> v.y >> v.z;
            vertices.emplace_back(v);
            aster_data >> l;
            ver_data << v.x << " " << v.y << " " << v.z << endl;
        }
        ver_data.close();
        cout << "number of vertices " << vertices.size() << endl;
        if (vertices.size() == 0)
            existence = false;
        while (l == 'f') {
            int i, j, k;
            aster_data >> i >> j >> k;
            polygons.emplace_back(polygon(vertices[i-1], vertices[j-1], vertices[k-1]));
            aster_data >> l;
        }
        cout << "number of polygons " << polygons.size() << endl;
        aster_data.close();
        for (auto& v:vertices) {
            if (v.dist() > maximum)
                maximum = v.dist();
            if (v.rms() > scale_factor)
                scale_factor = v.rms();
        }
    }

    bool is_inside(double x, double y, double z){
        bool k = false;
        for (auto& p:polygons)
            if (p.right_side(x, y, z))
                k = !(k);
        return k;
    }

//---------------- old --------------------------
//    bool is_inside(double x, double y, double z){
//        for (auto& p:polygons)
//            if (!p.right_side(x, y, z))
//                return false;
//        return true;
//    }

//--------------> uncomment <-------------------
    double step(){
        double min_norm = numeric_limits<double>::max();
        for (auto& p:polygons)
            if (p.area() < min_norm)
                min_norm = p.area();
        return 0.5*sqrt(2/sqrt(3)*min_norm);
    }

};

double side(double a, double R){
    int i=4, j=2;
    double s = numeric_limits<double>::max();
    double r_mer, r_par, r_pol;
    while (s >= a){
        j=(i+1)/2;
        r_mer = PI * R / (i-j+1);
//        if ((i-j+1)%2 == 0)
            r_par = PI / (2*j) * R * cos(PI / (i-j+1));
//        else
//            r_par = PI / (2*j) * R * cos(PI / 2 / (i-j+1));
        r_pol = PI / (2*j) * R * sin(PI / (i-j+1));
        s = r_par* r_mer;
        s = sqrt(s);
        i++;
    }
    return i;
}

template<typename T>
std::string toString(T val)
{
    std::ostringstream oss;
    oss << val;
    return string(oss.str());
};

template<typename T>
T fromString(const std::string& s)
{
  std::istringstream iss(s);
  T res;
  iss >> res;
  return res;
}

double dot_product(const state_vector &u, const state_vector &v){
    double result=0;
    for (int i=0;i<u.size();i++)
        result += u[i]*v[i];
    return result;
};

state_vector cross_product(const state_vector &u, const state_vector &v){
    state_vector result(u.size());
    result[0] = u[1]*v[2] - u[2]*v[1];
    result[1] = u[2]*v[0] - u[0]*v[2];
    result[2] = u[0]*v[1] - u[1]*v[0];
    return result;
};

double kron_delta(double n, double k){
    if (n == k)
        return 1;
    else
        return 0;
};

double sgn(double n){
    if (n < 0)
        return -1;
    else
        if (n > 0)
            return 1;
        else
            return 0;
};

state_matrix mult_mat(const state_matrix& A, const state_matrix& B){
    state_matrix C(3);
    for (auto& str:C)
        str.resize(3);
    for (int i=0; i<3;i++)
        for (int j=0;j<3;j++)
            for (int k=0;k<3;k++)
                C[i][j] += A[i][k]*B[k][j];
    return C;
};

state_matrix transpose_mat(const state_matrix& A){
    state_matrix C(A[0].size());
    for (auto& str:C)
        str.resize(A.size());
    for (int i=0; i<A[0].size();i++)
        for (int j=0;j<A.size();j++)
            C[i][j] = A[j][i];
    return C;
};

state_matrix matrix_rot(int elem_i, int elem_j, int dim, double theta){
    state_matrix A(dim);
    for (auto& str:A)
        str.resize(dim);
    for (int i=0;i<dim;i++)
        A[i][i] = 1;
    A[elem_i][elem_i] = cos(theta);
    A[elem_j][elem_j] = cos(theta);
    A[elem_j][elem_i] = sin(theta);
    A[elem_i][elem_j] = -sin(theta);
    return A;
};

class tensor{
    state_matrix matrix;
//public:
    void rot_m(state_matrix& m){
        m.resize(matrix.size());
        for (auto& mat:m)
            mat.resize(matrix.size());
        state_matrix m_new(m.size());
        for (auto& mat:m_new)
            mat.resize(m.size());
        for (int i=0;i<matrix.size();i++){
            for (int j=0;j<matrix.size();j++){
                if (i == j)
                    m[i][j] = 1;
                else
                    m[i][j] = 0;
     //           cout << m[i][j];
            }
   //         cout << endl;
        }
    //    cout << endl;
        double eps = numeric_limits<double>::epsilon();
        double t_m = 1.0, t_diag = 1.0;
        double theta, i_max, j_max, max_el;
        int count_it = 0;
        while (t_m/t_diag > eps){
            max_el=0;
            t_m = t_diag = 0;
            for (int i=0;i<matrix.size();i++){
                for (int j=i+1;j<matrix.size();j++){
                    if (fabs(matrix[i][j]) > max_el){
                        max_el = fabs(matrix[i][j]);
                        i_max = i;
                        j_max = j;
                    }
                    t_m += 2*matrix[i][j]*matrix[i][j];
                }
                t_diag += matrix[i][i]*matrix[i][i];
            }
            theta = atan2(2*matrix[i_max][j_max], matrix[i_max][i_max] - matrix[j_max][j_max])/2;
            m_new = matrix_rot(i_max, j_max, m.size(), theta);
            m = mult_mat(m, m_new);
            matrix = mult_mat(transpose_mat(m_new), matrix);
            matrix = mult_mat(matrix,m_new);
//            cout << "ok" << endl;
//            cout << "theta = " << theta << " i " << i_max << " j " << j_max << endl << endl;
//
//            cout << "m" << endl;
//            for (int i=0;i<3;i++){
//                for (int j=0;j<3;j++)
//                    cout << m_new[i][j] << " ";
//                cout << endl;
//            }
//            cout << endl;
//
//            state_vector u(3),v(3);
//            for (int j=0;j<2;j++){
//                for (int i=0;i<3;i++){
//                    u[i] = m_new[i][j];
//                    v[i] = m_new[i][j+1];
//                }
//                cout << dot_product(u,v) << endl;
//            }
//            cout << endl;
//
//            cout << "full m" << endl;
//            for (int i=0;i<3;i++){
//                for (int j=0;j<3;j++)
//                    cout << m[i][j] << " ";
//                cout << endl;
//            }
//            cout << endl;
//
//            for (int j=0;j<2;j++){
//                for (int i=0;i<3;i++){
//                    u[i] = m[i][j];
//                    v[i] = m[i][j+1];
//                }
//                cout << dot_product(u,v) << endl;
//            }
//            cout << endl;
//
//            cout << "tensor" << endl;
//            for (int i=0;i<3;i++){
//                for (int j=0;j<3;j++)
//                    cout << matrix[i][j] << " ";
//   //             cout << endl;
//            }
  //          cout << endl;
  //          count_it++;
        }

        //cout << count_it << endl;

        m = transpose_mat(m);

        std::ofstream vect_in_new, vect_in;
        vect_in_new.open("vect_in_new.dat",std::ios_base::out);
        for (int i=0;i<3;i++){
            vect_in_new << 0 << " " << 0 << " " << 0;
            for (int j=0;j<3;j++)
                vect_in_new << " " << matrix[i][j]/pow(10,13)*2500;
            vect_in_new << endl;
        }
        vect_in_new.close();

        vect_in.open("vect_in.dat",std::ios_base::out);
        for (int i=0;i<3;i++){
            vect_in << 0 << " " << 0 << " " << 0;
            for (int j=0;j<3;j++)
                vect_in << " " << m[i][j]*300;
            vect_in << endl;
        }
        vect_in.close();
    }
public:
//    tensor():matrix{{1,-2,3},{-2,1,-2},{3,-2,1}}{}
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
        std::ofstream vect_in_old;
        vect_in_old.open("vect_in_old.dat",std::ios_base::out);
        for (auto& ten:matrix){
            vect_in_old << 0 << " " << 0 << " " << 0;
            for (auto& elem:ten){
                elem *= mps.size();
                vect_in_old << " " << elem/pow(10,13)*2500;
            }
            vect_in_old << endl;
        }
        vect_in_old.close();
    }

    void transform(vector<m_point>& mps){
        state_matrix m(3);
        for(auto& ten:m)
            ten.resize(3);
        rot_m(m);
        for (auto& mp:mps){
            m_point copy_mp(mp);
            for (int i=0;i<3;i++){
                mp.coord[i] = 0;
                for (int j=0;j<3;j++)
                    mp.coord[i] += m[i][j] * copy_mp.coord[j];
            }
        }
    }
};

double Leg_pol(double n, double k, double phi) {
    double cos_phi=cos(phi), sin_phi=sin(phi);
    double pol_prev=0, pol_next, pol=1;

    for (int j=0;j<k;j++) {
        pol *= (2*j+1)*fabs(sin_phi)*sqrt((1+kron_delta(0,j))*static_cast<double>(2*j+3)/(2*j+1)/(2*j+1)/(2*j+2));
    }
    for (int i=k;i<n;i++) {
        pol_next = ((2*i+1)*cos_phi*pol*sqrt(static_cast<double>(2*i+3)/(2*i+1)*(i-k+1)/(i+k+1)) - (pol_prev == 0 ? 0:(i+k)*pol_prev*sqrt(static_cast<double>(2*i+3)/(2*i-1)*(i-k+1)/(i+k+1)*(i-k)/(i+k))))/(i-k+1);
        pol_prev = pol;
        pol = pol_next;
    }
    return pol;
};

double inv_alpha(double n, double k){
    double mult=1;
    for (int i=n-k+1;i<=n+k;i++)
        mult*=i;
    return mult/(2-kron_delta(0,k));
};

class coef{
private:
    state_matrix A, B;
public:
    coef(vector<m_point>& mps, body& aster, int num): A{num+1},B{num+1} {
        for (int i=0;i<=num;i++) {
            A[i].resize(i+1);
            B[i].resize(i+1);
            for (int j=0;j<=i;j++){
                A[i][j] = 0;
                B[i][j] = 0;
                for (auto& mp : mps){
                    A[i][j] += pow(mp.r(),i)*Leg_pol(i,j,mp.theta())*cos(j*mp.lambda());
                    B[i][j] += pow(mp.r(),i)*Leg_pol(i,j,mp.theta())*sin(j*mp.lambda());
                }
                A[i][j] /= inv_alpha(i,j)*pow(aster.scale_factor,i)*mps.size();
                B[i][j] /= inv_alpha(i,j)*pow(aster.scale_factor,i)*mps.size();
//                cout << i << " " << j << " " << A[i][j] << " " << B[i][j] << endl;
            }
//            cout << endl;
        }
    }

    void out_coef(const string& str, int num){
        ofstream file_coef(str);
        for (int i=0; i<=num; i++){
            for (int j=0; j<=i; j++)
                file_coef << i << " " << j << " " << A[i][j] << " " << B[i][j] << endl;
            file_coef << endl;
        }
        file_coef.close();
    }

    double iter_incr(double r, double tau, double f){
        return r + f*tau;
    };

    state_vector pol_coef(int theta, int lambda, int num){
        state_vector pc(num+1);
        for (int i=0; i<=num; i++)
            for (int j=0; j<=i; j++)
                pc[i] += Leg_pol(i,j,(double)theta*PI/180)*(A[i][j]*cos(j*(double)lambda*PI/180)+B[i][j]*sin(j*(double)lambda*PI/180));
        return pc;
    }

    double pol(int theta, int lambda, double r, double value, int num){
        state_vector pc = pol_coef(theta, lambda, num);
        double p = -value * pow(r,num+1);
//        cout << r << " " << p << endl;
        for (int i=0;i<=num;i++){
            p += pc[i] * pow(r,num-i);
//            cout << r << " " << p << endl;
        }
        return p;
    }

    double pol_der(int theta, int lambda, double r, double value, int num){
        state_vector pc = pol_coef(theta, lambda, num-1);
        double pd = -value*(num+1) * pow(r,num);
//        cout << r << " " << pd << endl;
        for (int i=0;i<=num-1;i++){
            pd += pc[i]*(num-i) * pow(r,num-i-1);
//            cout << r << " " << pd << endl;
        }
        return pd;
    }

    void dist(const string& str, double pot_value, int num_coef, int h_t, int h_l){
        double pol_coef, pol_v=1, tau;
        double surf_dist, surf_dist_prev=0;
        double x, y, z;
        double eps = numeric_limits<double>::epsilon()*10000;
        ofstream file_coef(str);
            for (int theta=0;theta<=180;theta+=h_t){
                for (int lambda=0;lambda<=360;lambda+=h_l){
                    surf_dist = 1/pot_value;
                    if (pol(theta, lambda, surf_dist, pot_value, num_coef) > 0){
                        while (pol(theta, lambda, surf_dist, pot_value, num_coef) > 0)
                            surf_dist += 0.01;
//                        cout << endl;
                    }
                    else {
                        while (pol(theta, lambda, surf_dist, pot_value, num_coef) < 0)
                            surf_dist -= 0.01;
//                        cout << endl;
                    }

                    tau = 1/fabs(pol_der(theta, lambda, surf_dist, pot_value, num_coef));
//                    cout << tau << endl;

//    //                    int pol_counter=0;
                    while (fabs(surf_dist - surf_dist_prev) > eps){
                        surf_dist_prev = surf_dist;
                        surf_dist = iter_incr(surf_dist,tau,pol(theta, lambda, surf_dist, pot_value, num_coef));
//                        cout << surf_dist << " " << pol(theta, lambda, surf_dist, pot_value, num_coef) << endl;
//    //                        if (pol_counter < 3)
//    //                            pol_counter++;
//    //                        else
//    //                            exit(0);
                    }
//                    cout << theta << " " << lambda << " " << surf_dist << endl;
                    x = surf_dist*sin(theta*PI/180)*cos(lambda*PI/180);
                    y = surf_dist*sin(theta*PI/180)*sin(lambda*PI/180);
                    z = surf_dist*cos(theta*PI/180);
                    file_coef << lambda << " " << theta << " " << surf_dist << " " << x << " " << y << " " << z <<  endl;

                }
                cout << "ok" << endl;
            }
        file_coef.close();
    }

};

int main(){
    double x, y, z, h, sum_pot=0, r, pot_value;
    int n, num_coef, ast_count;

    std::ifstream file_in, aster_data;

    file_in.open("file_in.dat",std::ios_base::in);
//  ввод параметров эллипса и шага
    file_in >> n >> num_coef >> ast_count >> pot_value;
    file_in.close();
    cout << n << " " << num_coef << " " << ast_count << endl;

    for (int counter = 100; counter <= 126; counter++){
        string str_c = toString(counter);
        cout << "asteroid " << counter << endl;

        body aster{};
        aster.open("./ast_"+str_c+"/"+str_c+".dat", "./ast_"+str_c+"/vertices.dat");
//        cout << "asteroid "+str_c << " " << aster.is_open() << endl;
        if (!aster.existence)
            cout << "file doesn't exist" << endl << endl;
        else {
            double m = aster.maximum;
            h = 2.0 * m / (n-1);
            //h = aster.step();
            cout << m << " " << h << endl;

//----------------------------------------------------------------------------------------------------------------------
            ofstream
                out_old_coord("./ast_"+str_c+"/"+str_c+"_old_coord_new.dat",std::ios_base::out),
                out_new_coord("./ast_"+str_c+"/"+str_c+"_new_coord_new.dat",std::ios_base::out);

            vector<m_point> disc_body{};
            disc_body.reserve(8*n*n*n);

            for (x=-m; x<=m; x+=h) {
                for (y=-m; y<=m; y+=h) {
                    for (z=-m; z<=m; z+=h) {
                        if (aster.is_inside(x, y, z)) {
        //                    out << x << " " << y << " " << z << endl;
                            disc_body.emplace_back(m_point(x,y,z));
                            for (int i=0; i<3; i++){
                                out_old_coord << m_point(x,y,z).coord[i] << " ";
                            }
                            out_old_coord << endl;
                        }
                    }
                }
            }

//            for (auto& p:disc_body){
//                for (int i=0;i<3;i++)
//                    out_old_coord << p.coord[i] << " ";
//                out_old_coord << endl;
//            }
            out_old_coord.close();
            cout << disc_body.size() << endl;

//----------------------------------------------------------------------------------------------------------------------
            state_vector center_m(3);
            for (auto& p:disc_body){
                for (int i=0;i<3;i++){
                    center_m[i] += p.coord[i]/disc_body.size();
                }
            }
        //    cout << center_m[0] << " " << center_m[1] << " " << center_m[2] << endl;
            for (auto& p:disc_body){
                for (int i=0;i<3;i++){
                    p.coord[i] -= center_m[i];
        //            out << p.coord[i] << " ";
                }
        //        out << endl;
            }

            tensor t{disc_body};
            t.transform(disc_body);

            for (auto& p:disc_body){
                for (int i=0;i<3;i++)
                    out_new_coord << p.coord[i] << " ";
                out_new_coord << endl;
            }
            out_new_coord.close();

            cout << "ok" << endl;
            //cout << side(h, m) << endl;

//----------------------------------------------------------------------------------------------------------------------
            coef st_coef(disc_body, aster, num_coef);
            st_coef.out_coef("./ast_"+str_c+"/"+str_c+"_coef.dat", num_coef);
            //st_coef.dist("./ast_"+str_c+"/"+str_c+"_eqvip.dat", pot_value, num_coef, 1, 3);
            cout << endl << "ok" << endl;

//----------------------------------------------------------------------------------------------------------------------
//            ofstream
//            res("./ast_"+str_c+"/"+str_c+"_res.dat",std::ios_base::out),
//            graph("./ast_"+str_c+"/"+str_c+"_graph.dat",std::ios_base::out);
//            r = 1.1*m;
//            state_vector prob_p(3);
//            prob_p[0] = 0; prob_p[1] = 0; prob_p[2] = -r;
//            for (int i=0;i<disc_body.size();i++) {
//                sum_pot += disc_body[i].calc_pot(prob_p);
//            }
//            res << "theta= " << 180 << " " << "lambda= " << 0 << " sum_pot= " << sum_pot << endl;
//            for (int theta=179;theta>0;theta-=1) {
//                for (int lambda=0; lambda<=360; lambda+=1) {
//                    sum_pot=0;
//                    prob_p[0]=r*sin((double) theta*PI/180.0)*cos((double) lambda*PI/180.0);
//                    prob_p[1]=r*sin((double) theta*PI/180.0)*sin((double) lambda*PI/180.0);
//                    prob_p[2]=r*cos((double) theta*PI/180.0);
//                    for (int i=0;i<disc_body.size();i++) {
//                        sum_pot += disc_body[i].calc_pot(prob_p);
//                    }
//                    res << " lambda= " << lambda << "theta= " << theta << " sum_pot= " << sum_pot << endl;
//                    graph << lambda << " " << theta << " " << sum_pot << endl;
//                }
//            }
//            sum_pot=0;
//            prob_p[0] = 0; prob_p[1] = 0; prob_p[2] = r;
//            for (int i=0;i<disc_body.size();i++) {
//                sum_pot += disc_body[i].calc_pot(prob_p);
//            }
//            res << "theta= " << 0 << " lambda= " << 0 << " sum_pot= " << sum_pot << endl;
//            res.close();
//            graph.close();
//----------------------------------------------------------------------------------------------------------------------
        }
    }

    return 0;
}


