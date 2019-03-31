//
//  SceneObject.h
//  informatiqueGraphique
//
//  Created by amr miguil on 24/03/2019.
//  Copyright © 2019 amr miguil. All rights reserved.
//

class SceneObject{
public:
    Vector centre;
    double rayon;
    Vector albedo;
    Vector partieSpeculaire;
    double n;
    Vector emissivite;
    TypeMateriel matType;
    Vector A;
    Vector B;
    Vector C;
    char obj;
    Vector offset;
    Vector normal;
    double redimension;
    
    SceneObject(const Vector centre , double rayon, Vector albedo, double n, TypeMateriel type = Diffus) : centre(centre), rayon(rayon), albedo(albedo), partieSpeculaire(partieSpeculaire), n(n), emissivite(emissivite), matType(type) {};
    SceneObject(const Vector A, const Vector B, const Vector C, const Vector albedo, double n, TypeMateriel type): A(A), B(B), C(C), albedo(albedo), partieSpeculaire(partieSpeculaire), n(n), matType(type)  {};
    SceneObject(const char obj, double redimension, Vector offset, Vector albedo, double n, TypeMateriel type): obj(obj), albedo(albedo), n(n), matType(type), redimension(redimension), offset(offset) {};
    virtual  bool isIntersect( Rayon &rayon, Vector &n, double &t, Vector &alphaBetaGamma, Vector &P, Vector &textures)    = 0;
};
