//
//  Triangle.h
//  informatiqueGraphique
//
//  Created by amr miguil on 04/03/2019.
//  Copyright © 2019 amr miguil. All rights reserved.





// matrice avec a b en première ligne puis c d
double determinant(const double &a, const double &b, const double &c, const double &d){
    return a * d - b * c;
}
class Triangle: public SceneObject{
    
public:
    Triangle(const Vector A, const Vector B, const Vector C, const Vector albedo, const Vector partieSpeculaire, const double n, const TypeMateriel type = Diffus): SceneObject(A,B,C,albedo, n, type) {
        
        
    };
    
    bool  isIntersect( Rayon &rayon, Vector &n, double &t, Vector &alphaBetaGamma, Vector &P, Vector &textures) {
        n = vectorProduct(B - A,  C - A);
        n.normalize();

        t = dot(C - rayon.origine, n) / dot(rayon.direction, n);
        if(t< 0) return false;
        P = rayon.origine + t * rayon.direction;
        
        double delta = determinant((B - A).norm2(), dot(C-A, B-A), dot(B-A, C-A), (C-A).norm2());
        double beta = determinant(dot(P-A, B - A), dot(C-A, B-A), dot(P-A, C-A), (C-A).norm2()) / delta;
        double gamma = determinant((B - A).norm2(), dot(P-A, B - A), dot(B-A, C-A),  dot(P-A, C-A)) / delta;
        
        if(beta >= 0 && beta <= 1 && gamma <= 1 && gamma >= 0 && (1 - beta - gamma) >= 0 && (1 - beta - gamma) <= 1){
            alphaBetaGamma = Vector((1 - beta - gamma), beta, gamma);
            return true;
        }
        
        return false;
    }
    
};

