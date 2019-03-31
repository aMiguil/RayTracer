//
//  Geometry.h
//  informatiqueGraphique
//
//  Created by amr miguil on 24/03/2019.
//  Copyright © 2019 amr miguil. All rights reserved.
//



class BoiteEnglobante{
public:
    BoiteEnglobante(){};
    BoiteEnglobante(const double xmin, const double xmax, const double ymin, const double ymax, const double zmin, const double zmax): xmin(xmin), ymin(ymin), xmax(xmax), ymax(ymax), zmin(zmin), zmax(zmax){};
    
    double xmin, xmax, ymin, ymax, zmin, zmax;
    // matrice a b première ligne puis c d
    double determinant(const double a, const double b, const double c, const double d){
        return a * d - b * c;
    }
    bool isIntersect(const Rayon &rayon) const{
        
        double t1x = (xmin - rayon.origine.x) / rayon.direction.x;
        double t2x = (xmax - rayon.origine.x) / rayon.direction.x;
        double t1y = (ymin - rayon.origine.y) / rayon.direction.y;
        double t2y = (ymax - rayon.origine.y) / rayon.direction.y;
        double t1z = (zmin - rayon.origine.z) / rayon.direction.z;
        double t2z = (zmax - rayon.origine.z) / rayon.direction.z;
        
        double txmin = std::min(t1x, t2x);
        double tymin = std::min(t1y, t2y);
        double tzmin = std::min(t1z, t2z);
        
        double txmax = std::max(t1x, t2x);
        double tymax = std::max(t1y, t2y);
        double tzmax = std::max(t1z, t2z);
        
        if(std::min(std::min(txmax, tymax), tzmax) < 0) return false;
        if(std::min(std::min(txmax, tymax), tzmax) - std::max(std::max(txmin, tymin), tzmin) > 0) return true;
        return false;
    }
    
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


class HierarchieBoitesEnglobantes{
    
public:
    int indiceInf, indiceSup;
    BoiteEnglobante be;
    HierarchieBoitesEnglobantes *gauche, *droit;
};

#include <map>
class Geometry : public SceneObject {
private:
    BoiteEnglobante be;
public:
    HierarchieBoitesEnglobantes bvh;
    
    Geometry(const char *obj, double redimension, Vector offset, Vector albedo, double n, TypeMateriel type): SceneObject(*obj, redimension, offset, albedo, n,type)  {
        readOBJ(obj);
        
        for (int i = 0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] * redimension + offset;
        }
        
        if(useBVH){
            construireHierarchie(&bvh, 0, indices.size());
        }else{
            be = construireBoiteEnglobante(0, vertices.size());
        }
        std::vector<std::string> texts;
        texts.push_back("visage.bmp");
        texts.push_back("cheveux.bmp");
        texts.push_back("corps.bmp");
        texts.push_back("pantalon.bmp");
        texts.push_back("accessoires.bmp");
        texts.push_back("mains.bmp");
        
        
        std::string prefix = pathToTexturesFolder;
        for(std::string path : texts){
            loadTextures(prefix + path);
        }
        std::cout<< "end resizing"<<std::endl;
    };
    
    void loadTextures(std::string prefix){
        
        const char *filename = prefix.c_str();
        textures.resize(textures.size() + 1);
        texturesWidths.resize(texturesWidths.size() + 1);
        texturesHeights.resize(texturesHeights.size() + 1);

        FILE* f;
        f = fopen(filename, "rb");
        unsigned char info[54];
        fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

        texturesWidths[texturesWidths.size() - 1] = *(int*)&info[18]; // extract image height and width from header
        texturesHeights[texturesHeights.size() - 1] = *(int*)&info[22];

        int size = 3 * texturesWidths[texturesWidths.size() - 1] * texturesHeights[texturesHeights.size() - 1];
        textures[textures.size() - 1].resize(size); // allocate 3 bytes per pixel
        fread(&textures[textures.size() - 1][0], sizeof(unsigned char), size, f); // read the rest of the data at once
        fclose(f);

        for (int i = 0; i < size; i += 3) {
            std::swap(textures[textures.size() - 1][i], textures[textures.size() - 1][i + 2]);
        }
    }
    
    BoiteEnglobante construireBoiteEnglobante(int indiceInf, int indiceSup){
        BoiteEnglobante temp;
        
        temp.xmin = vertices[indices[indiceInf].vtxi].x;
        temp.xmax = vertices[indices[indiceInf].vtxi].x;
        temp.ymin = vertices[indices[indiceInf].vtxi].y;
        temp.ymax = vertices[indices[indiceInf].vtxi].y;
        temp.zmin = vertices[indices[indiceInf].vtxi].z;
        temp.zmax = vertices[indices[indiceInf].vtxi].z;
        
        
        for(int i = indiceInf; i < indiceSup; i ++){
            temp.xmin = std::min(temp.xmin, vertices[indices[i].vtxi].x);
            temp.ymin = std::min(temp.ymin, vertices[indices[i].vtxj].y);
            temp.zmin = std::min(temp.zmin, vertices[indices[i].vtxk].z);
            temp.xmin = std::min(temp.xmin, vertices[indices[i].vtxj].x);
            temp.ymin = std::min(temp.ymin, vertices[indices[i].vtxk].y);
            temp.zmin = std::min(temp.zmin, vertices[indices[i].vtxi].z);
            temp.xmin = std::min(temp.xmin, vertices[indices[i].vtxk].x);
            temp.ymin = std::min(temp.ymin, vertices[indices[i].vtxi].y);
            temp.zmin = std::min(temp.zmin, vertices[indices[i].vtxj].z);
            
            temp.xmax = std::max(temp.xmax, vertices[indices[i].vtxi].x);
            temp.ymax = std::max(temp.ymax, vertices[indices[i].vtxj].y);
            temp.zmax = std::max(temp.zmax, vertices[indices[i].vtxk].z);
            temp.xmax = std::max(temp.xmax, vertices[indices[i].vtxj].x);
            temp.ymax = std::max(temp.ymax, vertices[indices[i].vtxk].y);
            temp.zmax = std::max(temp.zmax, vertices[indices[i].vtxi].z);
            temp.xmax = std::max(temp.xmax, vertices[indices[i].vtxk].x);
            temp.ymax = std::max(temp.ymax, vertices[indices[i].vtxi].y);
            temp.zmax = std::max(temp.zmax, vertices[indices[i].vtxj].z);
            
        }
        
        return temp;
    }
    void construireHierarchie(HierarchieBoitesEnglobantes *noeud, const int indiceInf, const int indiceSup){
        noeud->be = construireBoiteEnglobante(indiceInf, indiceSup);
        noeud->indiceInf = indiceInf;
        noeud->indiceSup = indiceSup;

        Vector beMin = Vector(noeud->be.xmin, noeud->be.ymin, noeud->be.zmin);
        Vector beMax = Vector(noeud->be.xmax, noeud->be.ymax, noeud->be.zmax);
        
        Vector diagonale = beMax - beMin;
        int dimensionLaPlusGrande;
        diagonale.x = std::abs(diagonale.x);
        diagonale.y = std::abs(diagonale.y);
        diagonale.z = std::abs(diagonale.z);
        if(diagonale.x > diagonale.y && diagonale.x > diagonale.z){
            dimensionLaPlusGrande = 0;
        }else if(diagonale.y > diagonale.x && diagonale.y > diagonale.z){
            dimensionLaPlusGrande = 1;
        }else{
            dimensionLaPlusGrande = 2;
        }
        /** On découpe la boite englobante selon la direction de la plus grande diagonale */
        double valeurDecoupage = beMin.getByInd(dimensionLaPlusGrande) + diagonale.getByInd(dimensionLaPlusGrande) * 0.5;
        int pivot = indiceInf-1;
        /** Boucle for pour parcourir les triangles qui se trouvent dans la boite englobante et on utilise l'algorithme quicksort pour trier les triangles selon droite ou gauche. La valeur du pivot finale sera la limite entre les triangles de droite et ceux de gauche  */
        for (int i =indiceInf; i < indiceSup; i++){
            /** On raisonne avec l'isobarycentre du triangle d'indice i. Et on prend toujours sa dimension donnée par dimensionLaPlusGrande  */
            double centerSplitDim = (vertices[indices[i].vtxi].getByInd(dimensionLaPlusGrande) + vertices[indices[i].vtxj].getByInd(dimensionLaPlusGrande) + vertices[indices[i].vtxk].getByInd(dimensionLaPlusGrande)) /3;
            
            /** Si ce point se trouve dans la demi-boite englobante donnée par le plan formé de dimensionLaPlusGrande = valeurDecoupage alors on déplace la valeur du pivot et on change les emplacements des triangles i et pivot. */
            if (centerSplitDim <= valeurDecoupage) {
                pivot++;

                std::swap(indices[i], indices[pivot]);
            }
        }
        /** Comme condition d'arrêt, une boite englobante à avec deux triangles. */
        if (pivot < indiceInf || pivot >= indiceSup- 1 || indiceInf + 1 == indiceSup){
            return;
        }
        
        noeud->gauche = new HierarchieBoitesEnglobantes();
        construireHierarchie(noeud->gauche, indiceInf, pivot+1);
        
        noeud->droit = new HierarchieBoitesEnglobantes();
        construireHierarchie(noeud->droit, pivot+1, indiceSup);
        
    }
    
    void readOBJ(const char* obj) {
        
        //char matfile[255];
        char grp[255];
        std::cout << obj << std::endl;
        FILE *f;
        f = fopen(obj, "r");
        int curGroup = -1;
        if (f==NULL) {
            printf("Error %d \n", errno);
        }
        std::map<std::string, int> groupNames;
        
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
                Vector vec;
                
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec.x, &vec.z, &vec.y, &col.x, &col.y, &col.z) == 6) {
                    col.x = std::min(1., std::max(0., col.x));
                    col.y = std::min(1., std::max(0., col.y));
                    col.z = std::min(1., std::max(0., col.z));
                    
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
                    
                } else {
                    //                    sscanf(line, "v %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
                    sscanf(line, "v %lf %lf %lf\n", &vec.x, &vec.z, &vec.y);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                //                sscanf(line, "vn %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
                sscanf(line, "vn %lf %lf %lf\n", &vec.x, &vec.z, &vec.y);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec.x, &vec.y);
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
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            //                                                        nn = sscanf(consumedline, "%u//%u %u//%u %u//%u\n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
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
                    } else {
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
                        } else {
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
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
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
    
//    std::vector<unsigned char*> textures;
    std::vector<std::vector<unsigned char> > textures;
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    std::vector<int> texturesWidths;
    std::vector<int> texturesHeights;
    
    
    bool isIntersect( Rayon &rayon, Vector &n, double &t, Vector &alphaBetaGamma, Vector &P, Vector &intersectedTextures)  {

        bool isIntersect = false;
        t = std::numeric_limits<double>::max();
        
        if(useBVH){
            if(!bvh.be.isIntersect(rayon)){
                return false;
            }
            
            std::list<const HierarchieBoitesEnglobantes*> boitesATraiter;
            boitesATraiter.push_back(&bvh);
            
            
            
            while (!boitesATraiter.empty()) {
                
                const HierarchieBoitesEnglobantes* current = boitesATraiter.back();
                boitesATraiter.pop_back();
                if (current->gauche && current->gauche->be.isIntersect(rayon)) {
                    
                        boitesATraiter.push_back(current->gauche);
                        
                }
                if (current->droit && current->droit->be.isIntersect(rayon)) {
                    
                        boitesATraiter.push_back(current->droit);
                        
                }
                
                if (!current->droit && !current->gauche) {
                    //
                    for(int i = current->indiceInf; i < current->indiceSup; i ++){
                        
                        Vector i0 = vertices[indices[i].vtxi];
                        Vector i1 = vertices[indices[i].vtxj];
                        Vector i2 = vertices[indices[i].vtxk];
                        
                        Triangle triangle = Triangle(i0, i1, i2, albedo, n, matType);
                        double tIntersTriangle;
                        Vector tempN, tempP, texts;
                        bool isIntersectTriangle = triangle.isIntersect(rayon, tempN, tIntersTriangle, alphaBetaGamma, tempP, texts);
                        if(isIntersectTriangle){
                            isIntersect = true;
                            if(tIntersTriangle < t){
                                t = tIntersTriangle;
                                n = normals[indices[i].ni] * alphaBetaGamma.x + normals[indices[i].nj] * alphaBetaGamma.y + normals[indices[i].nk] * alphaBetaGamma.z;
                                n.normalize();
                                P = tempP;
                                
                                int x = (uvs[indices[i].uvi].x * alphaBetaGamma.x + uvs[indices[i].uvj].x * alphaBetaGamma.y + uvs[indices[i].uvk].x * alphaBetaGamma.z) * (texturesWidths[indices[i].group] - 1);
                                int y = (uvs[indices[i].uvi].y * alphaBetaGamma.x + uvs[indices[i].uvj].y * alphaBetaGamma.y + uvs[indices[i].uvk].y * alphaBetaGamma.z) * (texturesHeights[indices[i].group] - 1);
                                
                                double cr = (textures[indices[i].group][(y*texturesWidths[indices[i].group] + x) * 3]) / 255.;
                                double cg = (textures[indices[i].group][(y*texturesWidths[indices[i].group] + x) * 3 + 1]) / 255.;
                                double cb = (textures[indices[i].group][(y*texturesWidths[indices[i].group]  + x) * 3 + 2]) / 255.;
                                
                                intersectedTextures = Vector(cr, cg, cb);
                            }
                            
                            
                        }
                        
                        
                    }
                }
                
            }
        }else{
            if(!be.isIntersect(rayon)){
                return false;
            }
            for(int i = 0; i < indices.size(); i +=1){
                
                Vector i0 = vertices[indices[i].vtxi];
                Vector i1 = vertices[indices[i].vtxj];
                Vector i2 = vertices[indices[i].vtxk];
                
                Triangle triangle = Triangle(i0, i1, i2, albedo, n, matType);
                double tIntersTriangle;
                Vector tempN, tempP, texts;
                bool isIntersectTriangle = triangle.isIntersect(rayon, tempN, tIntersTriangle, alphaBetaGamma, tempP, texts);
                if(isIntersectTriangle){
                    isIntersect = true;
                    if(tIntersTriangle < t){
                        t = tIntersTriangle;
                        n = normals[indices[i].ni] * alphaBetaGamma.x + normals[indices[i].nj] * alphaBetaGamma.y + normals[indices[i].nk] * alphaBetaGamma.z;
                        n.normalize();
                        P = tempP;
                        
                        int x = (uvs[indices[i].uvi].x * alphaBetaGamma.x + uvs[indices[i].uvj].x * alphaBetaGamma.y + uvs[indices[i].uvk].x * alphaBetaGamma.z) * (texturesWidths[indices[i].group] - 1);
                        int y = (uvs[indices[i].uvi].y * alphaBetaGamma.x + uvs[indices[i].uvj].y * alphaBetaGamma.y + uvs[indices[i].uvk].y * alphaBetaGamma.z) * (texturesHeights[indices[i].group] - 1);
                        
                        double red = (textures[indices[i].group][(y*texturesWidths[indices[i].group] + x) * 3]) / 255.;
                        double green = (textures[indices[i].group][(y*texturesWidths[indices[i].group] + x) * 3 + 1]) / 255.;
                        double blue = (textures[indices[i].group][(y*texturesWidths[indices[i].group]  + x) * 3 + 2]) / 255.;
                        
                        intersectedTextures = Vector(red, green, blue);
                    }
                    
                    
                }
                
                
            }
        }
        return isIntersect;
    }
    
    
    
};


