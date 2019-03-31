//
//  Sphere.h
//  informatiqueGraphique
//
//  Created by amr miguil on 24/03/2019.
//  Copyright © 2019 amr miguil. All rights reserved.
//



/**
 * intensityReflexion: how to reflect red, green and blue.
 */
class Sphere: public SceneObject{
    
public:
    Sphere( Vector centre , const double rayon, const Vector albedo, const double n, const TypeMateriel type = Diffus) :
    SceneObject(centre, rayon, albedo, n, type) {};
    
    bool isIntersect( Rayon &rayon, Vector &n, double &t, Vector &alphaBetaGamma, Vector &P, Vector &textures)   {
        Vector cMoinsO = rayon.origine -  this->centre ;
        double ac = cMoinsO.norm2() - this->rayon * this->rayon;
        
        t = std::numeric_limits<double>::max();
        const double delta = 4 * dot(rayon.direction, cMoinsO) * dot(rayon.direction, cMoinsO) - 4 * (ac);
        bool intersected = false;
        if (delta > 0) {
            double t1 = (-2 * dot(rayon.direction,cMoinsO) - sqrt(delta)) / 2;
            double t2 = (-2 * dot(rayon.direction,cMoinsO) + sqrt(delta)) / 2;
            
            if(t2 < 0){
                intersected =  false;
                return false;
            }if(t1 < 0){
                t = t2;
                intersected =  true;
                P = (rayon.origine + t * rayon.direction);
                n = P - centre;
                n.normalize();
            }
            if(t1 > 0){
                t = t1;
                intersected =  true;
                P = (rayon.origine + t * rayon.direction);
                n = P - centre;
                n.normalize();
            }
        } else if (delta == 0) {
            t =  -2 * dot(rayon.direction, rayon.origine - centre)/ 2;
            P = (rayon.origine + t * rayon.direction);
            n = P - centre;
            n.normalize();
            intersected =  true;
        }
        
        
        return intersected;
    }
    
    
};


bool operator==(const Sphere& a, const Sphere& b){
    if(a.centre == b.centre && a.rayon == b.rayon){
        return true;
    }
    return false;
}

bool operator!=(const Sphere& a, const Sphere& b){
    return !operator==(a,b);
}

