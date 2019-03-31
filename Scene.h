

class Scene{
    
public:
    Scene() {};
    
    std::vector<SceneObject*> listObjects;
    
    /** Je n'utilise pas addObject pour reconnaitre le type d'élément que j'envoie, et pour éviter de mettre const dans la définition de listObjects qui aura une influence sur la signature de isIntersect de SceneObject */
    
    void addSphere( Sphere &s) { listObjects.push_back(&s); }
    void addTriangle( Triangle &s) { listObjects.push_back(&s); }
    void addGeometry( Geometry &s) { listObjects.push_back(&s); }

    /**
     * Routine d'intersection avec tous les objets de la scène. On fait une boucle for sur chaque objet et on compare les éventuelles distances renvoyées par chaque élément. Pour un objet de type Geometry, on renvoie le vecteur de couleurs RGB textures. Pour tous les autres, même si on injecte texturesColor dans la fonction, on laisse sa valeur à Vector(-1,-1,-1)
     */
    bool intersect( Rayon &rayon, Vector &P, Vector &n, int &indiceSphere, Vector &texturesColor)  {
        bool isIntersect = false;
        double t = std::numeric_limits<double>::max();
        Vector tempN, tempP;
        for (int i = 0; i < this->listObjects.size(); i ++) {
            double distanceToRay;
            Vector alphaBetaGamma;
            Vector textures = Vector(-1,-1,-1);
            bool sphereIntersects ;
           
            sphereIntersects = listObjects[i]->isIntersect(rayon, tempN, distanceToRay, alphaBetaGamma, tempP, textures);
            
            if (sphereIntersects) {
                if (t > distanceToRay) {
                    t = distanceToRay;
                    n = tempN;
                    P = tempP;
                    indiceSphere = i;
                    texturesColor = textures;
                    isIntersect = true;
                }
            }
        }
        
        if (!isIntersect) return false;
        return true;
    }
    
};
