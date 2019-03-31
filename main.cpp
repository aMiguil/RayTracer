//
//  main.cpp
//  informatiqueGraphique
//
//  Created by amr miguil on 04/02/2019.
//  Copyright © 2019 amr miguil. All rights reserved.
//
#define _CRT_SECURE_NO_WARNINGS 1
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION


 bool useBVH = true;
 bool antiAliasing = false;
 bool indirectColor = false;
 bool sourceDeLumiereSpherique = false;
 bool animation = false;
char* pathToTexturesFolder;
char* pathGirl;
char* pathOut;
char **argvGlobal;

#include "main.h"

//const int W = 128;
//const int H = 128;
//const int W = 1024;
//const int H = 1024;
const int W = 512;
const int H = 512;
//const int W = 2056;
//const int H = 2056;
Scene scene = Scene();
const Vector centreCamera = Vector(0, 0, 55);
static std::vector<unsigned char> image(W * H * 3, 0);
const double fov = 60 * M_PI / 180;
const double tanfov = tan(fov / 2);
int nombreRayonIndirect = 200;
int I = 4500000;


/** Lumière pour une source ponctuelle */
const Vector L = Vector(-10, 20, 40);


/**
 * Cette fonction prend en argument un vecteur n et renvoie autour de ce vecteur un vecteur (une direction) aléatoire.
 */
Vector randomCos(const Vector n){
    double r1 = distrib(engine);
    double r2 = distrib(engine);

    double x = std::cos(2 * M_PI * r1) * std::sqrt(1 - r2);
    double y = std::sin(2 * M_PI * r1) * std::sqrt(1 - r2);
    double z = std::sqrt(r2);
    /** Création d'un vecteur aléatoire en fonction du numéro du thread utilisé. (Problème de concurrence lors de l'utilisation de distrib(engine). => Chaque Thread crée un distrib(engine) en fonction d'un seul moteur de génération.
     */
    Vector vectorAleatoire = Vector(distrib(listEngines[omp_get_thread_num()]),distrib(listEngines[omp_get_thread_num()]),distrib(listEngines[omp_get_thread_num()]));
    
    Vector tangent1 = vectorProduct(n, vectorAleatoire);
    tangent1.normalize();
    Vector tangent2 = vectorProduct(vectorAleatoire, n);

    /** vecteur random dans le repère global de la scène */
    return x * tangent1 + y * tangent2 + z * n;
}
/**
 * Fonction clé pour notre sujet. Elle prend en argument un rayon un nombre de rebonds maximum et la scène qui n'est pas mise en variable globale.
 */
Vector getColor(Rayon &rayon, int &rebondNum)
{
    if(rebondNum == 0){
        return Vector(0,0,0);
    }
    Vector P;
    Vector n;
    Vector resultat;
    int intersectSphere;
    /** Définition de notre sphère emettrice de la lumière. */
//    const SceneObject* sphereLumiere = scene.listObjects[0];
    Vector texturesColor = Vector(-1,-1,-1);
    
    
    const bool isIntersect = scene.intersect(rayon, P, n, intersectSphere, texturesColor);
    if(isIntersect){
        if(scene.listObjects[intersectSphere]->matType == Emissif  && sourceDeLumiereSpherique){
            return scene.listObjects[intersectSphere]->albedo * 2550 ;     //intensité * la couleur blanche
        }
        
        if( scene.listObjects[intersectSphere]->matType == Reflectif){
            if (rebondNum > 0)
            {
                rayon.reflechir(P, n);
                rebondNum -= 1;
                resultat = resultat + getColor(rayon, rebondNum);
            }
        }else if(scene.listObjects[intersectSphere]->matType == Refractif){
            if(rebondNum > 0){
                rayon.refracter(P, n, scene.listObjects[intersectSphere]->n);
                rebondNum -= 1;
                resultat = resultat + getColor(rayon, rebondNum);
            }
            
        }else if(scene.listObjects[intersectSphere]->matType == Fresnel){
            if(rebondNum > 0){
                rayon.fresneliser(P, n, scene.listObjects[intersectSphere]->n);
                rebondNum -= 1;
                resultat = resultat + getColor(rayon, rebondNum);
            }
            
        }else{
            
            
            
            // source sphérique
            if(sourceDeLumiereSpherique){

                Vector newDirection = P - scene.listObjects[0]->centre;     // Vecteur qui va du point d'intersection vers la lumière
                newDirection.normalize();
                Vector directionAleatoireLightToSphere;
                if(indirectColor){   //Si pas de indirectColor on n'envoie qu'un seul rayon vers le centre de la sphère de lumière, sinon on envoie un rayon vers une direction aléatoire autour
                    directionAleatoireLightToSphere = randomCos(newDirection);
                }else{
                    directionAleatoireLightToSphere= newDirection;
                }
                Vector pointAleatoire = directionAleatoireLightToSphere * scene.listObjects[0]->rayon + scene.listObjects[0]->centre;  //Point surface sphère lumineuse selon la direction aléatoire donnée par randomCos, ou bien le point d'intersection si pas de indirectColor
                Vector wI = pointAleatoire - P;
                double distanceToSphereLum = wI.norm2();
                wI.normalize();
                pointAleatoire = pointAleatoire+ 0.00001 * wI;
                
                P = P + 0.0001 * n;
                 Rayon secondIntersection = Rayon(P, wI);
                ////
                Vector pBefore;
                Vector nBefore;
                Vector bis = Vector(-1,-1,-1);
                int intersectSecondSphere;
                
                bool isSecondIntersection = scene.intersect(secondIntersection, pBefore, nBefore, intersectSecondSphere, bis);

                double distancePPprime = (pBefore - P).norm2();

                if ( isSecondIntersection && distancePPprime < distanceToSphereLum * 0.999)
                {
                    resultat = Vector(0,0,0);
                }
                else{
                    if(texturesColor.x < 0){
                        resultat = scene.listObjects[intersectSphere]->albedo * I / (4 * M_PI * distanceToSphereLum) * std::max(0.0, dot(wI, n)) * dot(directionAleatoireLightToSphere, wI * -1) / dot(newDirection, directionAleatoireLightToSphere);
                    }else {
                        resultat = texturesColor * I / (4 * M_PI * distanceToSphereLum) * std::max(0.0, dot(wI, n)) * dot(directionAleatoireLightToSphere, wI * -1) / dot(newDirection, directionAleatoireLightToSphere);
                    }
                }

                
                
            }else /** Eclairage ponctuel, même s'il est déjà la limite de l'éclairage sphérique lorsque la rayon de la sphère émettrice tend vers 0, j'ai préféré laisser ce cas.*/
            {
                
                double dCarre = (P-L).norm2();
                Vector PL = L - P;
                PL.normalize();
                P = P + 0.00001 * PL;
                Rayon secondIntersection = Rayon(P, PL);
                Vector pBefore;
                Vector nBefore;
                Vector bis;
                int intersectSecondSphere;
                bool isSecondIntersection = scene.intersect(secondIntersection, pBefore, nBefore, intersectSecondSphere, bis);
                Vector pPrime = pBefore;
                double distancePPprime = (pPrime - P).norm2();
                if ( isSecondIntersection && distancePPprime < dCarre * 0.99)
                {
                    resultat = Vector(0,0,0);
                }
                else{
                    
                    if(texturesColor.x < 0 ){
                        resultat = scene.listObjects[intersectSphere]->albedo  * std::max(0.0, dot(PL, n)) * I / dCarre;
                    }else{
                        resultat = texturesColor  * std::max(0.0, dot(PL, n)) * I / dCarre;
                    }
                }
                
            }
            
            
            if(indirectColor){
                Vector randomDirection = randomCos(n);
                Vector initPoint = P + 0.0001*randomDirection;
                
                Rayon aleatoire = Rayon(initPoint, randomDirection);
                rebondNum -= 1;
                if(texturesColor.x < 0){
                    resultat = resultat + multiplyByComps(getColor(aleatoire, rebondNum),scene.listObjects[intersectSphere]->albedo);
                }else{
                    resultat = resultat + multiplyByComps(getColor(aleatoire, rebondNum),texturesColor);
                }
            }
        }
    }
    
    return resultat;
}


Vector generateRay(const int i,const int j) {
    
    double x = distrib(engine), y = distrib(engine), R=sqrt(-2*log(x));
    double dx =  R*cos(2*M_PI*y)* 0.5;
    double dy =  R*sin(2*M_PI*y)* 0.5;
    Vector v = Vector( j-W/2+0.5+dx, i-H/2+0.5 + dy, -H / (2 * tanfov));
    v.normalize();
    return v;
};


/** Tourner le vecteur d'un angle angleRotation autour de axis*/
Vector getRotation(int axis, Vector vecteur, double angleRotation){
    switch (axis) {
        case 0:
            return Vector(vecteur.x, vecteur.y * std::cos(angleRotation) - vecteur.z * std::sin(angleRotation), vecteur.y * std::sin(angleRotation) + vecteur.z * std::cos(angleRotation));
        case 1:
            return Vector(vecteur.x * std::cos(angleRotation) - vecteur.z * std::sin(angleRotation), vecteur.y, vecteur.x * std::sin(angleRotation) + vecteur.z * std::cos(angleRotation));
        case 2:
            return Vector(vecteur.x * std::cos(angleRotation) - vecteur.y * std::sin(angleRotation), vecteur.x * std::sin(angleRotation) + vecteur.y * std::cos(angleRotation), vecteur.z);

        default:
            return NULL;
            break;
    }}

/** Pour générer un certain nombre d'images lors d'une rotation */
void generateRotation(){
    int counter = 0;
    double pas = M_PI / 25;
    for(double ang = 0; ang <= 2*M_PI; ang += pas){
        scene.listObjects[0]->centre = getRotation(1, scene.listObjects[0]->centre, pas);
        Vector cameraBis = getRotation(1, centreCamera , ang);

#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < H; i++)
        {
            
            for (int j = 0; j < W; j++)
            {
                
                Vector color = Vector();
                if(indirectColor){
                    for (int ray = 0; ray < nombreRayonIndirect; ray++)
                    {
                        Rayon rayon;
                        Vector direction;
                        
                        if(antiAliasing){
                            direction = generateRay(i,j);
                            direction = getRotation(1, direction , ang);
                            direction.normalize();
                        }else{
                            direction = Vector(0.5 + j - W / 2, 0.5 + i - H / 2, -H / (2 * tanfov));
                            direction = getRotation(1, direction , ang);
                            direction.normalize();
                        }
                        rayon = Rayon(cameraBis, direction);
                        int numRays = 5;
                        color = color + getColor(rayon, numRays);
                    }
                    color = color * (1. / nombreRayonIndirect);
                }else{
                    Vector direction = Vector(0.5 + j - W / 2, 0.5 + i - H / 2, -H / (2 * tanfov));
                    direction = getRotation(1, direction , ang);
                    direction.normalize();
                    
                    Rayon rayon = Rayon(cameraBis, direction);
                    int numRays = 5;
                    
                    color =    getColor(rayon, numRays);
                    
                }
                
                image[((H - i - 1) * W + j) * 3 + 0] = std::min(255. ,  std::pow(255 * color.x, 1./2.2 ) );
                image[((H - i - 1) * W + j) * 3 + 1] = std::min(255. ,  std::pow(255 * color.y, 1./2.2 ) );
                image[((H - i - 1) * W + j) * 3 + 2] = std::min(255. ,  std::pow(255 * color.z, 1./2.2 ) );
                
            }
            
            
        }
//        std::string location = "/Users/amrmiguil/Desktop/ToBeGiffed/image" + std::to_string(counter) + ".jpg";
        std::string location = pathOut + std::to_string(counter) + ".jpg";

        stbi_write_png(location.c_str(), W, H, 3, &image[0], 0); // @suppress("Invalid arguments")
        if(indirectColor){
            std::cout << "ended writing the image: ";
            std::cout << counter << std::endl;
        }

        counter ++;
    }
}

void generateOnePic(){
#pragma omp parallel for schedule(dynamic, 5)
    for (int i = 0; i < H; i++)
    {
        
        for (int j = 0; j < W; j++)
        {
            
            Vector color = Vector();
            if(indirectColor){
                for (int ray = 0; ray < nombreRayonIndirect; ray++)
                {
                    Rayon rayon;
                    Vector direction;
                    
                    if(antiAliasing){
                        direction = generateRay(i,j);
                        direction.normalize();
                    }else{
                        direction = Vector(0.5 + j - W / 2, 0.5 + i - H / 2, -H / (2 * tanfov));
                        direction.normalize();
                    }
                    rayon = Rayon(centreCamera, direction);
                    int numRays = 4;
                    color = color + getColor(rayon, numRays);
                }
                color = color * (1. / nombreRayonIndirect);
            }else{
                
                //Rayon rayon;
                Vector direction = Vector(0.5 + j - W / 2, 0.5 + i - H / 2, -H / (2 * tanfov));
                direction.normalize();
                
                Rayon rayon = Rayon(centreCamera, direction);
                int numRays = 5;
                color =    getColor(rayon, numRays);
                
            }
            
            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255. ,  std::pow(255 * color.x, 1./2.2 ) );
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255. ,  std::pow(255 * color.y, 1./2.2 ) );
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255. ,  std::pow(255 * color.z, 1./2.2 ) );
        }
        
        
    }
        std::string location = pathOut;

    stbi_write_png(location.c_str(), W, H, 3, &image[0], 0); // @suppress("Invalid arguments")
    
}


int main(int argc,char **argv){
    pathOut = argv[1];
    pathGirl = argv[2];
    pathToTexturesFolder = argv[3];
   
    // calcul du temps d'exécution, on lance le chrono ici
    auto start = std::chrono::high_resolution_clock::now();

    if(sourceDeLumiereSpherique){
        I = I * 10;
    }
    for(int i = 0; i < 12; i++){
        listEngines.push_back(engine);
    }
    
    Sphere sphereLum = Sphere(Vector(-10, 20, 40), 1, Vector(1., 1.,1.),  14, Emissif);
//    Sphere sphereLum = Sphere(Vector(0, 30,  0), 2, Vector(1., 1.,1.), Vector(1, 1, 1),  14, Emissif);
    Sphere s1 = Sphere(Vector(0, -1000, 0), 970, Vector(0, 0, 1), 0.00001);
    Sphere s2 = Sphere(Vector(0, 0, -1000), 900, Vector(0, 0.5, 0.019), 99);
    Sphere s3 = Sphere(Vector(0, 0, 1000), 940, Vector(0.5, 0.1, 0.1),  99);
    Sphere s4 = Sphere(Vector(0, 1000, 0), 940, Vector(0.4, 0, 0), 99);
    Sphere s5 = Sphere(Vector(-1000, 0, 0), 940, Vector(1, 0.764, 0.33),  99);
    Sphere s6 = Sphere(Vector(1000, 0, 0), 940, Vector(0.5, 0.5, 0.5),  99);
    
    
    /******************************************
     *    Trois sphères au centre seulement:
     */
//    Sphere center = Sphere(Vector(0, 0, 0), 10, Vector(1., 1.,1.),  1.4, Reflectif);
    Sphere center = Sphere(Vector(0, 0, 0), 15, Vector(1., 1.,1.),  1.4, Reflectif);
    Sphere oreilleMickey1 = Sphere(Vector(15, 15, 15), 5, Vector(0.0, 0.0,0.0),   1.4, Refractif);
    Sphere oreilleMickey2 = Sphere(Vector(-15, 15, 15), 5, Vector(0.5, 0.5,0.5), 1.2, Reflectif);
    
    /**********************************************/
     
    
    Geometry girl = Geometry(pathGirl, 30, Vector(0,-30,0), Vector(1, 1,1), 1.5,Diffus);

    
    if(sourceDeLumiereSpherique){
        scene.addSphere(sphereLum) ;
    }
    scene.addSphere(s1 );
    scene.addSphere( s3);
    scene.addSphere( s2);
    scene.addSphere( s5);
    scene.addSphere( s6);
    scene.addSphere( s4);
    scene.addSphere(oreilleMickey1);
    scene.addSphere(oreilleMickey2);

    scene.addGeometry(girl);

    
    if(animation){
        generateRotation();
    }else{
        generateOnePic();
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> duree = end - start;
    std::cout << duree.count() << std::endl;
    return 0;
};

