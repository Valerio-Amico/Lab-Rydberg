#include "stdafx.h"
#include <stdio.h>
#include <iostream>

using namespace std;

void calcola_distanza();
void calcola_volume();

int main(){
    int status;
    float velocita, durata, distanza;
    status = 0; //menu
    while(status==0){
        cout << "0 -> menu\n1 -> calcola distanza\n2 -> calcola volume\n -1 -> exit\n";
        cin >> status;
        if(status==1){
            calcola_distanza();
            status = 0;
        }else if(status==2){
            calcola_volume();
            status = 0;
        }
    }
    return 0;
}

void calcola_distanza(){
    float dist, velocita, durata;
    printf("inserisci velocità\n-> ");
    cin >> velocita;
    printf("inserisci durata\n-> ");
    cin >> durata;
    dist = velocita * durata;
    cout << "la durata è -> ";
    cout << dist << " \n";
    return;
}

void calcola_volume(){
    float lato, vol;
    printf("inserisci lato\n-> ");
    cin >> lato;
    vol =  lato * lato * lato;
    cout << "il volume è -> ";
    cout << vol << " \n";
    return;
}