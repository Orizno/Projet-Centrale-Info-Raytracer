#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdio.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <random>
#include <list>

#define M_PI 3.1415

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0, 1);


class Vecteur { //Définition d'une classe pour les vecteurs
public:
    explicit Vecteur(double x = 0, double y = 0, double z = 0) {
        coordonnees[0] = x;
        coordonnees[1] = y;
        coordonnees[2] = z;
    };
    double operator[](int i) const { return coordonnees[i]; };
    double& operator[](int i) { return coordonnees[i]; };
    double NormeAuCarre() { //norme au carré
        return coordonnees[0] * coordonnees[0] + coordonnees[1] * coordonnees[1] + coordonnees[2] * coordonnees[2];
    }
    Vecteur normaliser() { //pour normaliser un vecteur
        double norme = sqrt(NormeAuCarre());
        return Vecteur(coordonnees[0] / norme, coordonnees[1] / norme, coordonnees[2] / norme);
    }
private:
    double coordonnees[3];
};

Vecteur operator+(const Vecteur& a, const Vecteur& b) { //somme de deux vecteurs
    return Vecteur(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vecteur operator-(const Vecteur& a, const Vecteur& b) { //différence de deux vecteurs
    return Vecteur(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vecteur operator-(const Vecteur& a) { //inverse d'un vecteur
    return Vecteur(-a[0], -a[1], -a[2]);
}

Vecteur operator*(double a, const Vecteur& b) { //multiplication d'un vecteur pour un double
    return Vecteur(a * b[0], a * b[1], a * b[2]);
}

Vecteur operator*(const Vecteur& a, double b) { //idem
    return Vecteur(a[0] * b, a[1] * b, a[2] * b);
}

Vecteur operator/(const Vecteur& a, double b) { //division d'un vecteur pour un double
    return Vecteur(a[0] / b, a[1] / b, a[2] / b);
}

Vecteur prodVect(const Vecteur& a, const Vecteur& b) { //produit vectoriel entre 2 vecteurs
    return Vecteur(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);

}

double prodScal(const Vecteur& a, const Vecteur& b) { //produit scalaire entre 2 vecteurs
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vecteur prodTermeaterme(const Vecteur& a, const Vecteur& b) { //produit terme à terme entre 2 vecteur
    return Vecteur(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}


Vecteur aleatoire(const Vecteur& N) {
    double u1 = uniform(engine);// nombre aléatoire entre 0 et 1
    double u2 = uniform(engine);
    double x = cos(2 * M_PI * u1) * sqrt(1 - u2);
    double y = sin(2 * M_PI * u1) * sqrt(1 - u2);
    double z = sqrt(u2);
    Vecteur T1;
    if (N[0] < N[1] && N[0] < N[2]) {
        T1 = Vecteur(0, N[2], -N[1]);
    }
    else {
        if (N[1] < N[2] && N[1] < N[0]) {
            T1 = Vecteur(N[2], 0, -N[0]);
        }
        else {
            T1 = Vecteur(N[1], -N[0], 0);
        }
    }

    T1 = T1.normaliser();
    Vecteur T2 = prodVect(N, T1);
    return z * N + x * T1 + y * T2;

}



class Rayon { //Classe pour les rayons (caractérisés par leur centre et leur direction)
public:
    Rayon(const Vecteur& C, const Vecteur& u) : C(C), u(u) {
    }
    Vecteur C, u;
};

class Objet { //Peut représenter une sphère ou un maillage
public:
    Objet() {};

    virtual bool intersection(const Rayon& r, Vecteur& P, Vecteur& Normale, double& t, Vecteur& couleur) = 0;
    Vecteur albedo;
    bool Miroir, Transparent;
};

class Sphere : public Objet { //classe pour les sphères (rajout d'un attribut pour savoir si c'est un miroir ou pas, et un autre pour la transparence)
public:
    Sphere(const Vecteur& O, double R, const Vecteur& albedo, bool Miroir = false, bool Transparent = false) : O(O), R(R) {
        this->albedo = albedo;
        this->Miroir = Miroir;
        this->Transparent = Transparent;
    };
    bool intersection(const Rayon& r, Vecteur& P, Vecteur& N, double& t, Vecteur& couleur) {
        // résout a*t² + b*t + c = 0
        double a = 1;
        double b = 2 * prodScal(r.u, r.C - O);
        double c = (r.C - O).NormeAuCarre() - R * R;
        double delta = b * b - 4 * a * c;
        if (delta < 0) return false; //pas de solution, donc pas d'intersection

        double racinedelta = sqrt(delta);
        double t2 = (-b + racinedelta) / (2 * a);
        if (t2 < 0) return false; //pas de solution positive, donc pas d'intersection

        double t1 = (-b - racinedelta) / (2 * a);
        if (t1 > 0) // t est la plus petite valeur positive entre t1 et t2
            t = t1;

        else
            t = t2;

        P = r.C + t * r.u; // Point d'intersection
        N = (P - O).normaliser(); //Normale au point d'intersection
        couleur = this->albedo;





        return true;

    };
    Vecteur O;
    double R;

};


class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};

class BoundingBox {
public:
    bool intersection(const Rayon& r) {//intersection rayon boite
        double x1 = (mini[0] - r.C[0]) / r.u[0];
        double x2 = (maxi[0] - r.C[0]) / r.u[0];
        double xmin = std::min(x1, x2);
        double xmax = std::max(x1, x2);

        double y1 = (mini[1] - r.C[1]) / r.u[1];
        double y2 = (maxi[1] - r.C[1]) / r.u[1];
        double ymin = std::min(y1, y2);
        double ymax = std::max(y1, y2);

        double z1 = (mini[2] - r.C[2]) / r.u[2];
        double z2 = (maxi[2] - r.C[2]) / r.u[2];
        double zmin = std::min(z1, z2);
        double zmax = std::max(z1, z2);

        double max = std::min(xmax, std::min(ymax, zmax));
        double min = std::max(xmin, std::max(ymin, zmin));
        if (max < 0) return false;
        return max > min;
    }
    Vecteur mini, maxi;
};

class Noeud {
public:
    Noeud* fg, * fd;
    BoundingBox b;
    int debut, fin;

};


class TriangleMesh : public Objet {
public:
    ~TriangleMesh() {}
    TriangleMesh(const Vecteur& albedo, bool Miroir = false, bool Transparent = false) {
        this->albedo = albedo;
        this->Miroir = Miroir;
        this->Transparent = Transparent;
        BVH = new Noeud;
    };

    BoundingBox constrBB(int debut, int fin) {
        bb.mini = Vecteur(1E9, 1E9, 1E9);
        bb.maxi = Vecteur(-1E9, -1E9, -1E9);
        for (int i = debut; i < fin; i++) {
            for (int j = 0; j < 3; j++) {
                bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxi][j]);
                bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxi][j]);
                bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxj][j]);
                bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxj][j]);
                bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxk][j]);
                bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxk][j]);
            }
        }
        return bb;
    }

    void constrBVH(Noeud* n, int debut, int fin) {

        n->debut = debut;
        n->fin = fin;
        n->b = constrBB(n->debut, n->fin); // on commence par construire la BB de l'ensemble des triangles du noeud
        Vecteur longueurs = n->b.maxi - n->b.mini;
        int dim;
        if (longueurs[0] >= longueurs[1] && longueurs[0] >= longueurs[2]) { // on prend la dimension correspondant à la longueur la plus longue de la BB
            dim = 0;
        }
        else {
            if (longueurs[1] >= longueurs[0] && longueurs[1] >= longueurs[2]) {
                dim = 1;
            }
            else {
                dim = 2;
            }
        }
        double milieu = (n->b.mini[dim] + n->b.maxi[dim]) / 2; // on calcule le milieu de cette dimension
        int pivot = n->debut; // le pivot est au début le premier indice
        for (int i = n->debut; i < n->fin; i++) { // pour chaque triangle de l'ensemble
            double barycentre = (vertices[indices[i].vtxi][dim] + vertices[indices[i].vtxj][dim] + vertices[indices[i].vtxk][dim]) / 3; // on calcule son barycentre selon la dimension prise
            if (barycentre < milieu) { // s'il est avant le milieu 
                std::swap(indices[i], indices[pivot]);// alors on déplace le triangle à la place du pivot (il sera donc parmi l'ensemble du premier enfant)
                pivot++; //le pivot est incrémenté
            }
        }


        if (pivot == debut || pivot == fin || (fin - debut < 5)) return; // Si tous les triangles sont du même côtés, ou s'il n'y a plus que 5 triangles dans le BVH, on arrête cette branche du BVH

        n->fg = new Noeud;
        n->fd = new Noeud;

        constrBVH(n->fg, n->debut, pivot); // sinon on continue, avec les 2 nouveaux enfants (les triangles avant le pivot, et ceux après)
        constrBVH(n->fd, pivot, n->fin);

    }

    void readOBJ(const char* obj) {

        char matfile[255];
        char grp[255];

        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ') {
                Vecteur vec;

                Vecteur col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcouleurs.push_back(col);

                }
                else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vecteur vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vecteur vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char* consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                }
                else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    }
                    else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        }
                        else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    }
                    else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        }
                        else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            }
                            else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                }
                                else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }

            }

        }
        fclose(f);

    }
    bool intersection(const Rayon& r, Vecteur& P, Vecteur& Normale, double& t, Vecteur& couleur) {
        if (!BVH->b.intersection(r)) return false; //Si pas d'intersection avec la BVH, pas d'intersection avec le maillage
        t = 1E9;
        bool intersecte = false;
        std::list<Noeud*> l;
        l.push_back(BVH);
        while (!l.empty()) {//Tant qu'on n'a pas testé tout ce qui peut potentiellement intersecter le rayon
            Noeud* c = l.front(); // on prend le premier noeud 
            l.pop_front(); // (et on le retire de la liste)
            if (c->fg) { // s'il a des enfants
                if (c->fg->b.intersection(r)) { //on teste l'intersection de son enfant gauche avec le rayon, et si c'est le cas, on le met dans l pour pouvoir tester juste après l'intersection de ses enfants
                    l.push_front(c->fg);
                }
                if (c->fd->b.intersection(r)) {
                    l.push_front(c->fd);
                }
            }
            else { //s'il n'a pas d'enfant, alors on fait une intersection rayon triangle sur tous les triangles du noeud
                for (int i = c->debut; i < c->fin; i++) {
                    const Vecteur& A = vertices[indices[i].vtxi]; //on recupère les sommets du triangle
                    const Vecteur& B = vertices[indices[i].vtxj];
                    const Vecteur& C = vertices[indices[i].vtxk];
                    Vecteur e1 = B - A; //vecteur AB
                    Vecteur e2 = C - A; //vecteur AC
                    Vecteur N = prodVect(e1, e2); //Normale
                    Vecteur AC = r.C - A;
                    Vecteur ACvectU = prodVect(AC, r.u);
                    double UscalN = prodScal(r.u, N);
                    double beta = -prodScal(e2, ACvectU) / UscalN;
                    double gamma = prodScal(e1, ACvectU) / UscalN;
                    double alpha = 1 - beta - gamma;
                    double tobjet = -prodScal(AC, N) / UscalN;
                    if (beta >= 0 && gamma >= 0 && alpha >= 0 && beta <= 1 && gamma <= 1 && alpha <= 1 && tobjet > 0) { //Dans le cas d'une intersection
                        intersecte = true;
                        if (tobjet < t) { //si l'intersection est plus proche que celles existantes, on la prend en compte
                            t = tobjet;
                            Normale = alpha * normals[indices[i].ni] + beta * normals[indices[i].nj] + gamma * normals[indices[i].nk]; //on prend la normale définie dans le fichier obj
                            Normale = Normale.normaliser();
                            P = r.C + t * r.u;
                            Vecteur UV = alpha * uvs[indices[i].uvi] + beta * uvs[indices[i].uvj] + gamma * uvs[indices[i].uvk];
                            int W = Wtex[0];
                            int H = Htex[0];
                            UV = prodTermeaterme(UV, Vecteur(W, H, 0));
                            int uvx = UV[0] + 0.5;
                            int uvy = UV[1] + 0.5;
                            uvx = uvx % W;
                            uvy = uvy % H;
                            if (uvx < 0) uvx += W;
                            if (uvy < 0) uvy += H;
                            uvy = H - uvy - 1;
                            couleur = Vecteur(std::pow(textures[0][(uvy * W + uvx) * 3] / 255., 2.2),
                                std::pow(textures[0][(uvy * W + uvx) * 3 + 1] / 255., 2.2),
                                std::pow(textures[0][(uvy * W + uvx) * 3 + 2] / 255., 2.2));















                        }

                    }
                }
            }

        }
        return intersecte;
    }

    void loadTexture(const char* filename) {
        int W, H, C;
        unsigned char* texture = stbi_load(filename, &W, &H, &C, 3);
        Wtex.push_back(W); // largeur image
        Htex.push_back(H); //hauteur image
        textures.push_back(texture);
    }

    std::vector<TriangleIndices> indices;
    std::vector<Vecteur> vertices;
    std::vector<Vecteur> normals;
    std::vector<Vecteur> uvs;
    std::vector<Vecteur> vertexcouleurs;
    std::vector<unsigned char*> textures;
    std::vector<int> Wtex, Htex;
    BoundingBox bb;
    Noeud* BVH;

};

class Scene { // classe de la scene, qui peut comporter plusieurs sphères
public:
    Scene() {};
    bool intersection(const Rayon& r, Vecteur& P, Vecteur& N, Vecteur& albedo, double& t, bool& miroir, bool& transp, int& id) {
        t = 1E9;
        bool intersecte = false;
        for (int i = 0; i < objets.size(); i++) {// pour chacun des objets de la scène
            Vecteur Pobjet, Nobjet, Coulobjet;
            double tobjet;
            if (objets[i]->intersection(r, Pobjet, Nobjet, tobjet, Coulobjet) && tobjet < t) { //s'il y a une intersection plus proche que celles existantes, alors on prend en compte celle-ci et pas les autres
                intersecte = true;
                t = tobjet;
                P = Pobjet;
                N = Nobjet;
                albedo = Coulobjet;
                miroir = objets[i]->Miroir;
                transp = objets[i]->Transparent;
                id = i;

            }

        }
        return intersecte;
    }

    Vecteur obtenirCouleur(const Rayon& r, int rebond, bool dernierDiffus) { //pour obtenir la couleur
        if (rebond > 5) return Vecteur(0., 0., 0.); // si on dépasse 5 rebonds, on renvoit la couleur noire
        else {
            double t;
            bool miroir, transp;
            Vecteur P, N, albedo, couleur;
            int id;
            if (intersection(r, P, N, albedo, t, miroir, transp, id)) {// en cas d'intersection avec un objet de la scène
                if (id == 0) { //1ere sphère correspond à la lumière
                    if (rebond == 0 || !dernierDiffus) {
                        double Rl = dynamic_cast<Sphere*>(objets[0])->R;
                        return Vecteur(I, I, I) / (4 * M_PI * M_PI * Rl * Rl);
                    }
                    else
                        return Vecteur(0, 0, 0);
                }

                if (miroir) { // si c'est un miroir
                    Vecteur Directionreflechie = r.u - 2 * prodScal(r.u, N) * N; // direction de la réflexion
                    Rayon Rayonreflechi(P + 0.00001 * N, Directionreflechie); //rayon réfléchi, partant du point d'intersection et allant dans la direction réfléchie
                    return obtenirCouleur(Rayonreflechi, rebond + 1, false);
                }
                else {
                    if (transp) { // si  il est transparent

                        double n1 = 1., n2 = 1.4; //indices optiques
                        Vecteur N2 = N;
                        if (prodScal(r.u, N) > 0) { //si on sort de la sphère, on inverse les n et la normale
                            std::swap(n1, n2);
                            N2 = -N;
                        }
                        double angle = 1 - n1 * n1 / (n2 * n2) * (1 - prodScal(r.u, N2) * prodScal(r.u, N2));
                        if (angle < 0) { // si l'angle est plus petit que 0, il y a reflexion
                            Vecteur Directionreflechie = r.u - 2 * prodScal(r.u, N) * N;// direction de la réflexion
                            Rayon Rayonreflechi(P + 0.00001 * N, Directionreflechie);//rayon réfléchi, partant du point d'intersection et allant dans la direction réfléchie
                            return obtenirCouleur(Rayonreflechi, rebond + 1, false);
                        }
                        Vecteur Tt = n1 / n2 * (r.u - prodScal(r.u, N2) * N2); //Composante tangentielle de la direction de réfraction
                        Vecteur Tn = -sqrt(angle) * N2; //Composante normale de la direction de réfraction
                        Vecteur Directionrefractee = Tt + Tn; //direction de réfraction
                        Rayon Rayonrefracte(P - 0.0001 * N2, Directionrefractee); //rayon réfracté, partant du point d'intersection et allant dans la direction réfractée
                        return obtenirCouleur(Rayonrefracte, rebond + 1, false);
                    }
                    else {
                        double Rl = dynamic_cast<Sphere*>(objets[0])->R;
                        Vecteur Ol = dynamic_cast<Sphere*>(objets[0])->O;
                        Vecteur w = aleatoire((P - L).normaliser()); // Vecteur qui a plus de chance d'être dirigé vers la source de lumière.
                        Vecteur xp = w * Rl + Ol;
                        Vecteur Pxp = xp - P;
                        double normePxp = sqrt(Pxp.NormeAuCarre());
                        Pxp = Pxp / normePxp;
                        Vecteur Pombre, Nombre, albedoombre;
                        double tombre;
                        bool miroirombre, transombre;
                        int idombre;
                        Rayon Rayonombre(P + 0.00001 * N, Pxp);
                        if (intersection(Rayonombre, Pombre, Nombre, albedoombre, tombre, miroirombre, transombre, idombre) && tombre < normePxp - 0.0001) { //s'il y a intersection avant d'arriver à la source de lumière
                            couleur = Vecteur(0., 0., 0.);//pas éclairé = ombre
                        }
                        else {
                            couleur = I / (4 * M_PI * M_PI * Rl * Rl) * albedo / M_PI * std::max(0., prodScal(N, Pxp)) * std::max(0., -prodScal(w, Pxp)) / (normePxp * normePxp) / (std::max(0., -prodScal((L - P).normaliser(), w)) / (M_PI * Rl * Rl));
                        }

                        //eclaraige ind
                        Vecteur wi = aleatoire(N);
                        Rayon Rayonwi(P + 0.00001 * N, wi);//rayon aléatoire
                        couleur = couleur + prodTermeaterme(albedo, obtenirCouleur(Rayonwi, rebond + 1, true));


                    }
                }
            }
            return couleur;
        }
    }
    std::vector<Objet*> objets;
    Vecteur L;
    double I;

};


int main() {
    int W = 512;
    int H = 512;

    Vecteur C(0, 0, 55);
    Scene scene;
    scene.I = 5E9; //intensité lumineuse
    scene.L = Vecteur(-10, 20, 40);
    Sphere Slumiere(scene.L, 5, Vecteur(1., 1, 1));
    Sphere S1(Vecteur(0, 0, 0), 10, Vecteur(1., 0.3, 0.2));
    Sphere Sgauche(Vecteur(0, 0, 1000), 940, Vecteur(0., 1., 0.));
    Sphere Sdroite(Vecteur(0, 0, -1000), 940, Vecteur(1., 0., 1.));
    Sphere Shaut(Vecteur(0, 1000, 0), 940, Vecteur(1., 0., 0.));
    Sphere Sbas(Vecteur(0, -1000, 0), 990, Vecteur(0., 0., 1.));
    Sphere Splafond(Vecteur(1000, 0, 0), 940, Vecteur(1., 1., 1.));
    Sphere Ssol(Vecteur(-1000, 0, 0), 940, Vecteur(1., 1., 1.));
    TriangleMesh m(Vecteur(1., 1., 1.));
    m.readOBJ("Jeep_Renegade_2016.obj");
    m.loadTexture("car_jeep_ren.jpg");

    for (int i = 0; i < m.vertices.size(); i++) {
        m.vertices[i][0] *= 15;
        m.vertices[i][1] *= 15;
        m.vertices[i][2] *= 15;
        m.vertices[i][1] -= 10;
    }

    m.constrBVH(m.BVH, 0, m.indices.size());
    scene.objets.push_back(&Slumiere);
    scene.objets.push_back(&Sgauche);
    scene.objets.push_back(&Sdroite);
    scene.objets.push_back(&Shaut);
    scene.objets.push_back(&Sbas);
    scene.objets.push_back(&Splafond);
    scene.objets.push_back(&Ssol);
    scene.objets.push_back(&m);


    double fov = 60 * M_PI / 180;

    int nmbrayon = 50;

    std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for schedule(dynamic,1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vecteur couleur(0, 0, 0);
            for (int k = 0; k < nmbrayon; k++) {
                double u1 = uniform(engine);
                double u2 = uniform(engine);
                double x1 = 0.25 * cos(2 * M_PI * u1) * sqrt(-2 * log(u2));
                double x2 = 0.25 * sin(2 * M_PI * u1) * sqrt(-2 * log(u2));
                Vecteur u(j - W / 2 + x2 + 0.5, i - H / 2 + x1 + 0.5, -W / (2. * tan(fov / 2))); // antialiasing : on ne vise plus le centre des pixels, mais un endroit aléatoire proche de ce centre
                u = u.normaliser();

                double u3 = uniform(engine);
                double u4 = uniform(engine);
                double x3 = 0.01 * cos(2 * M_PI * u3) * sqrt(-2 * log(u4));
                double x4 = 0.01 * sin(2 * M_PI * u3) * sqrt(-2 * log(u4));

                Vecteur F = C + 55 * u; // point de focus
                Vecteur Cp = C + Vecteur(x3, x4, 0); // point d'origine aléatoire dans le disque d'ouverture
                Vecteur uversF = (F - Cp).normaliser(); // direction du rayon


                Rayon r(Cp, uversF);

                couleur = couleur + scene.obtenirCouleur(r, 0, false);
            }
            couleur = couleur / nmbrayon;


            couleur[0] = std::pow(couleur[0], 0.45);
            couleur[1] = std::pow(couleur[1], 0.45);
            couleur[2] = std::pow(couleur[2], 0.45);


            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., couleur[0]);
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., couleur[1]);
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., couleur[2]);
        }
    }
    stbi_write_png("image_textureJeep.png", W, H, 3, &image[0], 0);

    return 0;
}
