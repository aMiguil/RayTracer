//
//  classes.h
//  informatiqueGraphique
//
//  Created by amr miguil on 04/02/2019.
//  Copyright © 2019 amr miguil. All rights reserved.
//

/**
 * Type de matériel sous la forme d'une énumération
 */
enum TypeMateriel {
    Diffus, Reflectif, Refractif, Fresnel, Emissif
};


class Vector {
public:
    Vector(const double x = 0, const double y = 0, const double z = 0) :
    x(x), y(y), z(z) {
    }
    ;
    double x, y, z;
    
    /**
     * carré de la norme
     */
    double norm2() {
        return x * x + y * y + z * z;
    };

    void normalize() {
        double n = sqrt(norm2());
        x = x / n;
        y = y / n;
        z = z / n;
    }
    
    /**
     *   La fonction suivante sera utilisée lors de la construction de la hiérarchie des boites englobantes. Elle donne soit x soit y soit z en fonction de l'indice ind.
     */
    const double getByInd(const int ind){
        switch(ind){
            case 0:
                return this->x;
                break;
            case 1:
                return this->y;
                break;
            case 2:
                return this->z;
                break;
            default:
                break;
        }
        return -1;
    }
};
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}
;
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}
;
Vector operator*(const double &a, const Vector& b) {
    return Vector(a * b.x, a * b.y, a * b.z);
}
;
Vector operator*(const Vector& a, const double &b) {
    return Vector(a.x * b, a.y * b, a.z * b);
}
;
Vector operator/(const Vector& a, const double &b) {
    return Vector(a.x / b, a.y / b, a.z / b);
}
;

double dot(const Vector& a, const Vector& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

bool operator==(const Vector& a, const Vector& b) {
    if (a.x == b.x && a.y == b.y && a.z == b.z) {
        return true;
    }
    return false;
}

/**
 * Définition d'un produit de vecteurs composante par composante
 */
Vector multiplyByComps(const Vector &a, const Vector& b){
    return Vector(a.x * b.x, a.y * b.y, a.z * b.z);
}

/**
 * Produit vectoriel
 */
Vector vectorProduct(const Vector& a, const Vector& b){
    return Vector(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}




class Rayon {
public:
    Rayon(){};
    
    Rayon(const Vector origine , Vector direction) {
        this->origine = origine;
        this->direction = direction;
        //        direction.normalize();
    }
    ;
    Vector origine;
    Vector direction;
    /**
     * reflechir selon la formule de l'angle réfléchi sur une surface à un point P de normale n
     */
    void reflechir(const Vector P, const Vector n){
        this->direction = this->direction - 2 * dot(this->direction, n)* n;
        origine = P + 0.001*this->direction;
        
    }
    
    /**
     * réfracter selon la formule en faisant attention à la convention in/out et les réaffectations qui s'en découlent
     */
    void refracter(Vector P,  Vector n, double nSphere, double nAir = 1){
        double temp = -1;
        if (dot(direction, n) > 0){  // Si on sort de la sphère
            n = n* -1;
            temp = nSphere;
            nSphere = nAir;
            nAir = temp;
        }
        double racineRefract = 1 - (nAir / nSphere)*(nAir / nSphere)*(1 - dot(direction, n)*dot(direction, n));
        
        if(racineRefract >= 0){
            direction = (nAir / nSphere) * direction - ((nAir / nSphere) * dot(direction, n) + sqrt(racineRefract)) * n;
            origine = P + 0.001 * direction;
        }
        else{
            reflechir(P, n);
        }
    }
    
    /**
     * fonction dont le nom évoque le comportement à suivre dans le cas d'une surface avec effet "Fresnel".
     */
    void fresneliser(Vector P, Vector n, double nSphere, double nAir=1){
        double temp = -1;
        if (dot(direction, n) > 0){
            n = n* -1;
            temp = nSphere;
            nSphere = nAir;
            nAir = temp;
        }
        double k0 = std::pow((nSphere - nAir),2) / std::pow((nSphere + nAir),2);
        double r = k0 + (1-k0) * std::pow((1+ dot(n, direction)), 5);
        
        // On tire une valeur aléatoire entre 0 et 1 et en fonction du coefficient de réflexion on sera à une probabilité près en réflexion ou en réfraction. 
        double alea = distrib(engine);
        if(alea < r){
            reflechir(P,n);
        }else{
            
            double racineRefract = 1 - (nAir / nSphere)*(nAir / nSphere)*(1 - dot(direction, n)*dot(direction, n));
            
            if(racineRefract >= 0){
                direction = (nAir / nSphere) * direction - ((nAir / nSphere) * dot(direction, n) + sqrt(racineRefract)) * n;
                origine = P + 0.001 * direction;
            }
            else{
                reflechir(P, n);
            }
        }
    }

    
    
};




