#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#ifdef __unix__
#include <alloca.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#define HFLAG 0      // Remove Hydrogenes
#define SIMPLEFLAG 0 // Less is more
/*
 DockRMSD: an open-source tool for atom mapping and RMSD calculation of symmetric molecules through graph isomorphism

 Written by Eric Bell And reviewed by Kevin Barre
 
 v1.0 written 5/2/2019
 from latest update (v1.1) written 8/26/2019
 ######################################################################
 # DockRMSD (v1.1): docking pose distance calculation                 #
 #                                                                    #
 # If you use DockRMSD in your work, please cite:                     #
 #                                                                    #
 # Bell, E.W., Zhang, Y. DockRMSD: an open-source tool for atom       #
 # mapping and RMSD calculation of symmetric molecules through graph  #
 # isomorphism. Journal of Cheminformatics, 11:40 (2019).             #
 ######################################################################
*/

const int MAXBONDS = 6;        //Maximum number of bonds allowable on a single atom
const int MAXLINELENGTH = 150; //Maximum length (in characters) of a line in a mol2 file
const double MAXMAPCOUNT = 0;  //Maximum amount of possible mappings before symmetry heuristic is used
const int MAXDEPTH = 2;

int grabAtomCount(FILE *mol2, int hflag);
int arrayIdentity(char **arr1, char **arr2, int arrlen);
int inArray(int n, int *arr, int arrlen);
void readMol2(char **atoms, double **coords, char ***bonds, int *nums, FILE *mol2, int atomcount, int hflag);
int generalizeBonds(char ***bonds, int atomcount);
char **buildTree(int depth, int index, char **atoms, char ***bonds, char *prestring, int prevind, int atomcount);
double searchAssigns(int atomcount, int **allcands, int candcounts[], int *assign, char ***tempbond, char ***querybond, double **querycoord, double **tempcoord, int *bestassign);
struct DockRMSD assignAtoms(char **tempatom, char ***tempbond, char **queryatom, char ***querybond, double **querycoord, double **tempcoord, int *querynums, int *tempnums, int atomcount, int simpleflag);
int validateBonds(int *atomassign, int proposedatom, int assignpos, char ***querybond, char ***tempbond, int atomcount);
struct DockRMSD make_and_send_point(FILE *query, FILE *template);

typedef struct DockRMSD
{
    double rmsd;
    double total_of_possible_mappings;
    char *optimal_mapping;
    char *error;
} DockRMSD;

struct DockRMSD dock_rmsd(FILE *query, FILE *template)
{
    int querycount = grabAtomCount(query, HFLAG);
    int tempcount = grabAtomCount(template, HFLAG);
    struct DockRMSD rmsd = {0, 0, "", ""};
    if (querycount != tempcount)
    {
        rmsd.error = "Error: Query and template don't have the same atom count!";
        return rmsd;
    }
    if (querycount == 0)
    {
        rmsd.error = "Error: Query file has no atoms!";
        return rmsd;
    }
    if (tempcount == 0)
    {
        rmsd.error = "Error: Template file has no atoms!";
        return rmsd;
    }
    //Initialize pointer arrays and fill them with information from the mol2 files
    int i, j;
    char **queryatoms = alloca(sizeof(char *) * querycount);
    double **querycoords = alloca(sizeof(double *) * querycount);
    char ***querybonds = alloca(sizeof(char **) * querycount);
    char **tempatoms = alloca(sizeof(char *) * tempcount);
    double **tempcoords = alloca(sizeof(double *) * tempcount);
    char ***tempbonds = alloca(sizeof(char **) * tempcount);
    int *querynums = alloca(sizeof(int) * querycount);
    int *tempnums = alloca(sizeof(int) * tempcount);
    for (i = 0; i < querycount; i++)
    {
        char *queryatom = alloca(sizeof(char) * 3);
        *(queryatoms + i) = queryatom;
        char *tempatom = alloca(sizeof(char) * 3);
        *(tempatoms + i) = tempatom;
        double *querycoord = alloca(sizeof(double) * 3);
        *(querycoords + i) = querycoord;
        double *tempcoord = alloca(sizeof(double) * 3);
        *(tempcoords + i) = tempcoord;
        char **querybondrow = alloca(sizeof(char *) * querycount);
        char **tempbondrow = alloca(sizeof(char *) * tempcount);
        for (j = 0; j < querycount; j++)
        {
            char *querybond = alloca(sizeof(char) * 3);
            strcpy(querybond, "");
            *(querybondrow + j) = querybond;
            char *tempbond = alloca(sizeof(char) * 3);
            strcpy(tempbond, "");
            *(tempbondrow + j) = tempbond;
        }
        *(querybonds + i) = querybondrow;
        *(tempbonds + i) = tempbondrow;
    }
    readMol2(queryatoms, querycoords, querybonds, querynums, query, querycount, HFLAG);
    readMol2(tempatoms, tempcoords, tempbonds, tempnums, template, tempcount, HFLAG);
    fclose(query);
    fclose(template);
    if (!arrayIdentity(queryatoms, tempatoms, querycount))
    {
        rmsd.error = "Template and query don't have the same atoms.";
        return rmsd;
    }

    char **flatquerybonds = malloc(sizeof(char *) * (querycount * querycount));
    char **flattempbonds = malloc(sizeof(char *) * (tempcount * tempcount));
    for (i = 0; i < querycount; i++)
    {
        memcpy(flatquerybonds + querycount * i, *(querybonds + i), sizeof(char *) * querycount);
        memcpy(flattempbonds + tempcount * i, *(tempbonds + i), sizeof(char *) * tempcount);
    }
    if (!arrayIdentity(flatquerybonds, flattempbonds, querycount * querycount))
    {
        // Remove bond typing if they don't agree between query and template
        generalizeBonds(querybonds, querycount);
        generalizeBonds(tempbonds, tempcount);
        for (i = 0; i < querycount; i++)
        {
            memcpy(flatquerybonds + querycount * i, *(querybonds + i), sizeof(char *) * querycount);
            memcpy(flattempbonds + tempcount * i, *(tempbonds + i), sizeof(char *) * tempcount);
        }
        // If the general bonds still don't agree, the molecules aren't the same
        if (!arrayIdentity(flatquerybonds, flattempbonds, querycount * querycount))
        {
            rmsd.error = "Template and query don't have the same bonding network.";
            return rmsd;
        }
    }
    free(flatquerybonds);
    free(flattempbonds);
    return assignAtoms(tempatoms, tempbonds, queryatoms, querybonds, querycoords, tempcoords, querynums, tempnums, querycount, SIMPLEFLAG);
}

//Comparator for compatibility with qsort
int strcompar(const void *a, const void *b) { return strcmp(*(char **)a, *(char **)b); }

//Returns 1 if two arrays contain the same string elements, otherwise returns 0
int arrayIdentity(char **arr1, char **arr2, int arrlen)
{
    int i;
    char **list1 = (char **)alloca(sizeof(char *) * arrlen);
    char **list2 = (char **)alloca(sizeof(char *) * arrlen);
    for (i = 0; i < arrlen; i++)
    {
        list1[i] = *(arr1 + i);
        list2[i] = *(arr2 + i);
    }
    qsort(list1, arrlen, sizeof(list1[0]), strcompar);
    qsort(list2, arrlen, sizeof(list2[0]), strcompar);
    for (i = 0; i < arrlen; i++)
    {
        if (strcmp(list1[i], list2[i]))
        {
            return 0;
        }
    }
    return 1;
}

//Returns the index+1 if the element n is already in the array, otherwise returns 0
int inArray(int n, int *arr, int arrlen)
{
    int i;
    for (i = 0; i < arrlen; i++)
    {
        if (*(arr + i) == n)
        {
            return i + 1;
        }
    }
    return 0;
}

//Returns the count of atoms in a mol2 file
int grabAtomCount(FILE *mol2, int hflag)
{
    char line[MAXLINELENGTH];
    int atomcount = 0;
    int countflag = 0;
    while (fgets(line, MAXLINELENGTH, mol2) != NULL)
    {
        if (strlen(line) > 1 && line[strlen(line) - 2] == '\r')
        { //For windows line endings
            line[strlen(line) - 2] = '\n';
            line[strlen(line) - 1] = '\0';
        }
        if (!strcmp(line, "@<TRIPOS>ATOM\n"))
        {
            countflag = 1;
            continue;
        }
        if (!strcmp(line, "@<TRIPOS>BOND\n"))
        {
            countflag = 0;
            break;
        }
        if (countflag && strlen(line) > 1)
        {
            char *token = strtok(line, " \t");
            int i;
            for (i = 0; i < 5; i++)
            {
                token = strtok(NULL, " \t");
            }
            if (hflag || strcmp(token, "H"))
            {
                atomcount++;
            }
        }
    }
    if (ferror(mol2))
    {
        printf("Error %d while reading in file.\n", ferror(mol2));
    }
    rewind(mol2); //resets the file pointer for use in other functions
    return atomcount;
}

//Fills atoms, coords, and bonds with information contained within a mol2 file
void readMol2(char **atoms, double **coords, char ***bonds, int *nums, FILE *mol2, int atomcount, int hflag)
{
    int i = 0;
    int sectionflag = 0; //Value is 1 when reading atoms, 2 when reading bonds, 0 before atoms, >2 after bonds
    char line[MAXLINELENGTH];
    int *atomnums = (int *)alloca(sizeof(int) * atomcount);
    // int atomnums[atomcount]; //Keeps track of all non-H atom numbers for bond reading
    while (fgets(line, MAXLINELENGTH, mol2) != NULL)
    {
        if (strlen(line) > 1 && line[strlen(line) - 2] == '\r')
        { //Handling windows line endings
            line[strlen(line) - 2] = '\n';
            line[strlen(line) - 1] = '\0';
        }
        if (!strcmp(line, "@<TRIPOS>ATOM\n") || (sectionflag && line[0] == '@'))
        {
            sectionflag++;
        }
        else if (sectionflag == 1)
        { //Reading in atoms and coordinates
            double coord[3];
            int j = 0;
            char *parts = strtok(line, " \t");
            int atomnum = atoi(parts);
            parts = strtok(NULL, " \t");
            for (j = 0; j < 3; j++)
            {
                parts = strtok(NULL, " \t");
                coord[j] = atof(parts);
            }
            parts = strtok(NULL, " \t");
            if (hflag || strcmp("H", parts))
            {
                char *element = strtok(parts, ".");
                strcpy(*(atoms + i), element);
                atomnums[i] = atomnum;
                for (j = 0; j < 3; j++)
                {
                    *(*(coords + i) + j) = coord[j];
                }
                i++;
            }
        }
        else if (sectionflag == 2)
        { //Reading in bonding network
            char *parts = strtok(line, " \t");
            parts = strtok(NULL, " \t");
            int from = inArray(atoi(parts), atomnums, atomcount) - 1;
            parts = strtok(NULL, " \t");
            int to = inArray(atoi(parts), atomnums, atomcount) - 1;
            parts = strtok(NULL, " \t");
            parts = strtok(parts, "\n");
            if (from >= 0 && to >= 0)
            {
                strcpy(*(*(bonds + to) + from), parts);
                strcpy(*(*(bonds + from) + to), parts);
            }
        }
    }
    for (i = 0; i < atomcount; i++)
    {
        *(nums + i) = atomnums[i];
    }
}

//Changes all bond types to generic "b" if the bond types don't agree between query and template. Returns true if this has already been done, false if not.
int generalizeBonds(char ***bonds, int atomcount)
{
    int i;
    int j;
    for (i = 0; i < atomcount; i++)
    {
        for (j = 0; j < atomcount; j++)
        {
            if (strcmp(*(*(bonds + i) + j), ""))
            {
                if (!strcmp(*(*(bonds + i) + j), "b"))
                {
                    return 1;
                }
                strcpy(*(*(bonds + i) + j), "b");
            }
        }
    }
    return 0;
}

//Recursive function that returns a pointer array of leaves in the bonding tree at a specified depth
char **buildTree(int depth, int index,
                 char **atoms, char ***bonds,
                 char *prestring, int prevind, int atomcount)
{
    // int n;
    if (depth == 0)
    { //Base case, if max depth is reached return the prestring
        char **outp = (char **)malloc(sizeof(char *) * 2);
        char *outstring = (char *)malloc(sizeof(char) * (strlen(prestring) + 1));
        strcpy(outstring, prestring);
        *outp = outstring;
        *(outp + 1) = NULL;
        return outp;
    }
    else
    {
        //Grab all immediate neighbors of the current atom
        int bondinds[MAXBONDS];
        char *bondtypes[MAXBONDS];
        int i;
        int bondi = 0;
        for (i = 0; i < atomcount; i++)
        {
            if (strcmp(*(*(bonds + index) + i), ""))
            {
                bondinds[bondi] = i;
                bondtypes[bondi] = *(*(bonds + index) + i);
                bondi++;
            }
        }
        int maxleaves = (int)pow((double)MAXBONDS, (double)depth); //Maximum possible number of leaf nodes at this depth
        char **outlist = (char **)malloc(sizeof(char *) * (maxleaves + 1));
        if (!outlist)
        { //If the outlist pointer wasn't allocated, you've hit the recursion limit
            return NULL;
        }
        *outlist = NULL;
        int leafind = 0;
        for (i = 0; i < bondi; i++)
        {
            if (bondinds[i] != prevind)
            { //Don't analyze the atom we just came from in the parent function call
                char newpre[strlen(prestring) + 8];
                strcpy(newpre, prestring);
                strcat(newpre, bondtypes[i]);
                strcat(newpre, *(atoms + bondinds[i]));
                //Recurse and fetch all leaves of the binding tree for this neighbor
                char **new = buildTree(depth - 1, bondinds[i], atoms, bonds, newpre, index, atomcount);
                char **newit = new;
                while (*newit)
                { //Append the new leaves onto outlist
                    *(outlist + leafind) = *newit;
                    newit++;
                    leafind++;
                }
                free(new);
                *(outlist + leafind) = NULL;
            }
        }
        if (!*outlist)
        { //Another base case, if the current atom's only neighbor is the atom analyzed in the parent function call
            free(outlist);
            outlist = malloc(sizeof(char *) * 2);
            *outlist = malloc(sizeof(char) * (strlen(prestring) + 1));
            strcpy(*outlist, prestring);
            *(outlist + 1) = NULL;
        }
        return outlist;
    }
}

double searchAssigns(int atomcount, int **allcands,
                     int candcounts[], int *assign,
                     char ***tempbond, char ***querybond,
                     double **querycoord, double **tempcoord, int *bestassign)
{
    int i;
    int j;
    int index;
    double **dists = (double **)malloc(sizeof(double *) * atomcount); //Distances between query atoms and template atoms
    double *querydists = (double *)alloca(sizeof(double) * atomcount * atomcount);
    int **queryconnect = (int **)alloca(sizeof(int *) * atomcount * MAXBONDS);
    int *bondcount = (int *)alloca(sizeof(int) * atomcount);
    int *connectcount = (int *)alloca(sizeof(int) * atomcount);

    //precalculate all query-template atomic distances
    for (i = 0; i < atomcount; i++)
    {
        *(assign + i) = -1;
        connectcount[i] = 0;
        double *distind = (double *)malloc(sizeof(double) * candcounts[i]);
        for (j = 0; j < candcounts[i]; j++)
        {
            double dist = 0.0;
            for (index = 0; index < 3; index++)
            {
                dist += pow(*(*(querycoord + i) + index) - *(*(tempcoord + *(*(allcands + i) + j)) + index), 2.0);
            }
            *(distind + j) = dist;
        }
        *(dists + i) = distind;
    }
    memcpy(bestassign, assign, sizeof(int) * atomcount);

    //precalculate all query-query atomic distances for feasibility check
    for (i = 0; i < atomcount; i++)
    {
        for (j = i; j < atomcount; j++)
        {
            double dist = 0.0;
            for (index = 0; index < 3; index++)
            {
                dist += pow(*(*(querycoord + i) + index) - *(*(querycoord + j) + index), 2.0);
            }
            querydists[i * atomcount + j] = pow(dist, 0.5);
            querydists[j * atomcount + i] = querydists[i * atomcount + j];
        }
    }

    //Calculate bond degree for every atom
    for (i = 0; i < atomcount; i++)
    {
        int degree = 0;
        for (j = 0; j < atomcount; j++)
        {
            if (strcmp(*(*(querybond + i) + j), ""))
            {
                queryconnect[i][degree] = j;
                degree++;
            }
        }
        bondcount[i] = degree;
    }

    //bubble sort all possible atoms at each position by query-template distance
    for (index = 0; index < atomcount; index++)
    {
        for (i = 0; i < candcounts[index]; i++)
        {
            for (j = 0; j < candcounts[index] - i - 1; j++)
            {
                if (*(*(dists + index) + j) > *(*(dists + index) + j + 1))
                {
                    int tempind = *(*(allcands + index) + j);
                    double tempdist = *(*(dists + index) + j);
                    *(*(allcands + index) + j) = *(*(allcands + index) + j + 1);
                    *(*(dists + index) + j) = *(*(dists + index) + j + 1);
                    *(*(allcands + index) + j + 1) = tempind;
                    *(*(dists + index) + j + 1) = tempdist;
                }
            }
        }
    }
    int *history = (int *)alloca(sizeof(int) * atomcount);
    int *histinds = (int *)alloca(sizeof(int) * atomcount);
    index = 0;
    for (i = 0; i < atomcount; i++)
    {
        histinds[i] = 0;
    }

    double runningTotal = 0.0;
    double bestTotal = DBL_MAX;
    while (1)
    { //While not all mappings have been searched
        if (index == atomcount)
        { //If we've reached the end of a mapping and haven't been pruned
            if (runningTotal < bestTotal)
            {
                memcpy(bestassign, assign, sizeof(int) * atomcount);
                bestTotal = runningTotal;
            }
            index--;
            runningTotal -= *(*(dists + history[index]) + histinds[index] - 1);
            continue;
        }
        if (histinds[index])
        { //Atom to analyze has been picked, we need to change the mapping

            while (index > 0 && histinds[index] == candcounts[history[index]])
            {
                histinds[index] = 0;
                *(assign + history[index]) = -1;
                for (i = 0; i < bondcount[history[index]]; i++)
                {
                    connectcount[queryconnect[history[index]][i]]--;
                }
                index--;
                runningTotal -= *(*(dists + history[index]) + histinds[index] - 1);
            }
            if (index == 0 && histinds[0] == candcounts[history[0]])
            { //This occurs when all mappings have been exhausted
                break;
            }
        }
        else
        { //Pick an atom to analyze
            int nextAtom = 0;
            double bestMetric = DBL_MAX;
            for (i = 0; i < atomcount; i++)
            {
                double nodescore = -0.1 * ((double)connectcount[i]) + 1.0 * ((double)candcounts[i]);
                if (nodescore < bestMetric && *(assign + i) == -1)
                {
                    nextAtom = i;
                    bestMetric = nodescore;
                }
            }
            history[index] = nextAtom;
            for (i = 0; i < bondcount[history[index]]; i++)
            {
                connectcount[queryconnect[history[index]][i]]++;
            }
            //printf("Next atom: %d\n",nextAtom);
        }
        int foundflag = 0;
        for (i = histinds[index]; i < candcounts[history[index]]; i++)
        {

            if (runningTotal + *(*(dists + history[index]) + i) > bestTotal)
            { //Dead end elimination check
                break;
            }

            if (!inArray(*(*(allcands + history[index]) + i), assign, atomcount) && validateBonds(assign, *(*(allcands + history[index]) + i), history[index], querybond, tempbond, atomcount))
            { //Feasibility check
                foundflag = 1;
                *(assign + history[index]) = *(*(allcands + history[index]) + i);
                histinds[index] = i + 1;
                runningTotal += *(*(dists + history[index]) + i);
                index++;
                break;
            }
        }
        if (!foundflag)
        { //This occurs if none of the remaining possibilities can be mapped
            if (index == 0)
            {
                break;
            }
            else
            {
                histinds[index] = 0;
                *(assign + history[index]) = -1;
                for (i = 0; i < bondcount[history[index]]; i++)
                {
                    connectcount[queryconnect[history[index]][i]]--;
                }
                index--;
                runningTotal -= *(*(dists + history[index]) + histinds[index] - 1);
            }
        }
    }
    for (i = 0; i < atomcount; i++)
    {
        free(*(dists + i));
    }
    free(dists);
    if (*bestassign != -1)
    {
        return pow(bestTotal / ((double)atomcount), 0.5);
    }
    else
    {
        return DBL_MAX;
    }
}

//Checks if the assignment of the current atom is feasible
int validateBonds(int *atomassign, int proposedatom,
                  int assignpos, char ***querybond,
                  char ***tempbond, int atomcount)
{
    int j;
    int queri = 0;
    int queryinds[MAXBONDS];
    for (j = 0; j < atomcount; j++)
    {
        if (strcmp(*(*(querybond + assignpos) + j), ""))
        {
            queryinds[queri] = j;
            queri++;
        }
    }
    for (j = 0; j < queri; j++)
    {
        int assignatom = *(atomassign + queryinds[j]);
        if (assignatom >= 0 && !strcmp(*(*(tempbond + proposedatom) + assignatom), ""))
        {
            return 0;
        }
    }

    return 1;
}

//Returns the lowest RMSD of all possible mappings for query atoms with template indices given the two molecules' bonding network
DockRMSD assignAtoms(char **tempatom, char ***tempbond,
                     char **queryatom, char ***querybond,
                     double **querycoord, double **tempcoord,
                     int *querynums, int *tempnums, int atomcount,
                     int simpleflag)
{
    int i;
    struct DockRMSD rmsd = {0, 0, "", ""};
    int **allcands = (int **)malloc(sizeof(int *) * atomcount);  //List of all atoms in the template that could feasibly be each query atom
    int *candcounts = (int *)alloca(sizeof(char *) * atomcount); //Number of atoms in the template that could feasibly be each query atom

    //Iterate through each query atom and determine which template atoms correspond to the query
    for (i = 0; i < atomcount; i++)
    {
        int j;
        int candidates[atomcount]; //Flags corresponding to if each template atom could correspond to the current query atom
        int viablecands = 0;       //Count of template atoms that could correspond to the current query atom
        for (j = 0; j < atomcount; j++)
        {
            if (!strcmp(*(queryatom + i), *(tempatom + j)))
            {
                candidates[j] = 1;
                viablecands++;
            }
            else
            {
                candidates[j] = 0;
            }
        }
        int treedepth = 1; //Recursion depth
        while (treedepth <= MAXDEPTH)
        { //Recurse deeper until you've searched all atoms or you've hit the recursion limit
            char **qtree = buildTree(treedepth, i, queryatom, querybond, *(queryatom + i), -1, atomcount);
            if (!qtree)
            { //This means you've hit the recursion limit and can't analyze any deeper
                break;
            }
            char **qit = qtree;
            int qlen = 0;
            while (*qit)
            {
                qlen++;
                qit++;
            }
            for (j = 0; j < atomcount; j++)
            {
                if (candidates[j])
                {
                    char **ttree = buildTree(treedepth, j, tempatom, tempbond, *(tempatom + j), -1, atomcount);
                    char **tit = ttree;
                    int tlen = 0;
                    while (*tit)
                    {
                        tlen++;
                        tit++;
                    }
                    if (tlen != qlen || !arrayIdentity(qtree, ttree, qlen))
                    { //If the template atom tree and query atom tree don't have the same leaves, they're not the same atom
                        candidates[j] = 0;
                        viablecands--;
                    }
                    tit = ttree;
                    while (*tit)
                    {
                        free(*tit);
                        tit++;
                    }
                    free(ttree);
                }
            }
            qit = qtree;
            while (*qit)
            {
                free(*qit);
                qit++;
            }
            free(qtree);
            treedepth++;
        }
        if (!viablecands)
        { //If there's no possible atom, something went wrong or the two molecules are not identical
            if (!generalizeBonds(querybond, atomcount))
            {
                if (!simpleflag)
                {
                    char *formatstring = NULL;
                    asprintf(&formatstring, "No atoms mappable for atom %d, generalizing bonds...\n", i);
                    rmsd.error = formatstring;
                }
                generalizeBonds(tempbond, atomcount);
                for (j = 0; j < i; j++)
                {
                    free(*(allcands + j));
                    candcounts[j] = 0;
                }
                i = -1;
                continue;
            }
            else
            {
                char *formatstring = NULL;
                asprintf(&formatstring, "Atom assignment failed for atom %d.\n", i);
                rmsd.error = formatstring;
                return rmsd;
            }
        }
        else
        { //Otherwise, store all possible template atoms for this query atom
            candcounts[i] = viablecands;
            int *atomcands = malloc(sizeof(int) * viablecands);
            int k = 0;
            for (j = 0; j < atomcount; j++)
            {
                if (candidates[j])
                {
                    *(atomcands + k) = j;
                    k++;
                    if (k == viablecands)
                    {
                        break;
                    }
                }
            }
            *(allcands + i) = atomcands;
        }
    }

    double possiblemaps = 1.0;
    for (i = 0; i < atomcount; i++)
    {
        possiblemaps *= candcounts[i];
    }
    //Calculate RMSD of all possible mappings given each query atoms possible template atoms and return the minimum
    int *assign = (int *)malloc(sizeof(int) * atomcount);
    int *bestassign = (int *)malloc(sizeof(int) * atomcount);
    double bestrmsd = searchAssigns(atomcount, allcands, candcounts, assign, tempbond, querybond, querycoord, tempcoord, bestassign);
    char *header = "Optimal mapping (First file -> Second file, * indicates correspondence is not one-to-one):\n";
    char *optimal_mapping = (char *)malloc(strlen(header) + 1);
    char *formatstring = NULL;
    strcpy(optimal_mapping, header);
    rmsd.rmsd = bestrmsd;
    rmsd.total_of_possible_mappings = possiblemaps;
    if (bestrmsd == DBL_MAX)
    {
        rmsd.error = "No valid mapping exists\n";
        return rmsd;
    }
    for (i = 0; i < atomcount; i++)
    {
        if (0 > asprintf(&formatstring, "%s%3d -> %s%3d ", *(queryatom + i), *(querynums + i), *(tempatom + *(bestassign + i)), *(tempnums + *(bestassign + i))))
            return rmsd;
        optimal_mapping = (char *)realloc(optimal_mapping, strlen(optimal_mapping) + strlen(formatstring) + 1);
        strcat(optimal_mapping, formatstring);
        free(formatstring);
        if (*(querynums + i) == *(tempnums + *(bestassign + i)))
        {
            optimal_mapping = (char *)realloc(optimal_mapping, strlen(optimal_mapping) + 2);
            strcat(optimal_mapping, "\n");
        }
        else
        {
            optimal_mapping = (char *)realloc(optimal_mapping, strlen(optimal_mapping) + 3);
            strcat(optimal_mapping, "*\n");
        }
    }
    free(assign);
    free(bestassign);
    rmsd.optimal_mapping = optimal_mapping;
    return rmsd;
}
