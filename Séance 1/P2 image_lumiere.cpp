#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.1415

class Vecteur{ //D?finition d'une classe pour les vecteurs
public:
    explicit Vecteur(double x=0, double y=0, double z=0){
        coordonnees[0] = x;
        coordonnees[1] = y;
        coordonnees[2] = z;
    };
    double operator[](int i) const {return coordonnees[i]; };
    double &operator[](int i) {return coordonnees[i]; };
    double NormeAuCarre() { //norme au carr?
        return coordonnees[0]*coordonnees[0]+ coordonnees[1]*coordonnees[1]+ coordonnees[2]*coordonnees[2];
    }
    Vecteur normaliser() { //pour normaliser un vecteur
        double norme = sqrt(NormeAuCarre());
        return Vecteur(coordonnees[0] / norme,coordonnees[1] / norme,coordonnees[2] / norme);
    }
private:
    double coordonnees[3];
};

Vecteur operator+(const Vecteur& a, const Vecteur& b) { //somme de deux vecteurs
    return Vecteur(a[0] + b[0],a[1] + b[1],a[2] + b[2]);
}

Vecteur operator-(const Vecteur& a, const Vecteur& b) { //diff?rence de deux vecteurs
    return Vecteur(a[0] - b[0],a[1] - b[1],a[2] - b[2]);
}

Vecteur operator*(double a, const Vecteur& b) { //multiplication d'un vecteur pour un double
    return Vecteur(a*b[0],a*b[1],a*b[2]);
}

Vecteur operator*(const Vecteur& a, double b) { //idem
    return Vecteur(a[0]*b,a[1]*b,a[2]*b);
}

Vecteur operator/(const Vecteur& a, double b) { //division d'un vecteur pour un double
    return Vecteur(a[0]/b,a[1]/b,a[2]/b);
}


double prodScal(const Vecteur& a, const Vecteur& b) { //produit scalaire entre 2 vecteurs
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


class Rayon{ //Classe pour les rayons (caract?ris?s par leur centre et leur direction)
public:
    Rayon(const Vecteur& C, const Vecteur& u) : C(C), u(u){
    }
    Vecteur C,u;
};

class Sphere{ //classe pour les sph?res
public:
    Sphere(const Vecteur& O, double R): O(O), R(R) {
    }
    bool intersection(const Rayon& r, Vecteur& P, Vecteur& N){
        // r?sout a*t? + b*t + c = 0
        double a=1;
        double b=2*prodScal(r.u, r.C - O);
        double c= (r.C-O).NormeAuCarre() - R*R;
        double delta = b*b - 4 * a*c;
        if (delta < 0) return false; //pas de solution, donc pas d'intersection

        double racinedelta = sqrt(delta);
        double t2 = (-b + racinedelta)/(2*a);
        if (t2 < 0) return false; //pas de solution positive, donc pas d'intersection

        double t1 = (-b - racinedelta)/(2*a);
        double t; //plus petite valeur positive entre t1 et t2
        if (t1 > 0)
            t=t1;
        else
            t=t2;

        P = r.C + t*r.u; // Point d'intersection
        N = (P-O).normaliser(); //Normale au point d'intersection

        return true;

    };
    Vecteur O;
    double R;
};

int main() {
	int W = 512;
	int H = 512;

	Vecteur C(0, 0, 55);
    Vecteur O(0,0,0);
    double R = 10; //rayon de la sph?re
    Sphere S(O,R);
    double fov = 60 * M_PI / 180;
    double I = 1E7; //intensit? lumineuse
    Vecteur rho (1,0,0);
    Vecteur L(-10,20,40); //source lumineuse

	std::vector<unsigned char> image(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

            Vecteur u(j-W/2,i-H/2, -W / (2.*tan(fov / 2)));
            u=u.normaliser();
            Rayon r(C, u);
            Vecteur P,N ;
            Vecteur couleur(0,0,0);
            if (S.intersection(r,P,N)) { //en cas d'intersection avec la sph?re, on calcule la couleur du pixel
                double normePL = sqrt((L-P).NormeAuCarre()); //norme de PL
				double prov = std::max(0., prodScal(N,(L-P)/normePL));
                couleur = I/(4*M_PI*normePL*normePL) *rho/M_PI * prov;
            }

			image[((H - i - 1)*W + j) * 3 + 0] = std::min(255.,couleur[0]);
			image[((H - i - 1)*W + j) * 3 + 1] = std::min(255.,couleur[1]);
			image[((H - i - 1)*W + j) * 3 + 2] = std::min(255.,couleur[2]);
		}
	}
	stbi_write_png("image_lumi?re.png", W, H, 3, &image[0], 0);

	return 0;
}
