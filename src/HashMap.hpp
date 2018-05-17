/*
 * HashMap.hpp
 * 
 * Copyright (C) 2013 Louis-Francois Handfield
 * e-mail: lfhandfield@gmail.com
 *
 * This program is free software; upon notification by email to the licensor
 * of the licencee identity and nature of use, the licencee can redistribute
 * this program and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 
 * of the License, or (at the licencee option) any later version. As such,
 * no further notifications are required.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
	
#undef LFHTEMP
#define LFHTEMP template<class Key, class Data, class HashFunc >
LFHTEMP Data& myHashmap<Key,Data,HashFunc>::operator[](const Key & key){
    unsigned int j = find(key); 
    if (j == 0xFFFFFFFF){
        
        pair< KeyElem<Key, Data> , unsigned int> new_buck;
        new_buck.first.k = key;
        ExOp::toZero(new_buck.first.d);
       	int rsize = (heap.asize & 0x7FFFFFFF);
        j = heap.size();
        if ((((heap.asize + 1)^(heap.asize))> rsize)&&((heap.asize & 0x80000000) == 0)) { // going to up_alloc 
            heap.push_back(new_buck);
            rehash(hash_mag + 1);
        }else{
            heap.push_back(new_buck);
            unsigned int i = hashpos( HashFunc::makeSeed(key) );
            heap[j].second = hash[i];
            hash[i] = j;
        }

    }
    return heap[j].first.d; 
}

LFHTEMP Data myHashmap<Key,Data,HashFunc>::operator[](const Key & key) const{
    unsigned int j = find(key); 
    if (j == 0xFFFFFFFF) {Data fout; ExOp::toZero(fout); return fout;}
    return heap[j].first.d; 
}


LFHTEMP unsigned int myHashmap<Key,Data,HashFunc>::find(const Key & k )const{
    if (hash_mag == 0) return(0xFFFFFFFF);
    unsigned int j = hash[hashpos( HashFunc::makeSeed(k) )];
    for(;j != 0xFFFFFFFF;j = heap[j].second){
        if (k == heap[j].first.k) break;
    }
    return(j);
}


LFHTEMP void myHashmap<Key,Data,HashFunc>::swap_and_pop(unsigned int ite){
    // swap and pop
   // printf("pop swap %i, ,%i\n", ite, heap.size());fflush(stdout);
    int rsize = (heap.asize & 0x7FFFFFFF);
    if (ite == heap.size()-1){
    if ((((heap.asize)^(heap.asize-1))> rsize-1)&&(( (heap.asize & 0x80000000) != 0)||(rsize == 1) )  ) {
        heap.pop_back();rehash((rsize == 1) ? 0 : hash_mag-1);
    }else{
        heap.pop_back();
    }

    }else{

    heap[ite] = heap[heap.size()-1];
    
    
    
    
    
    if ((((heap.asize)^(heap.asize-1))> rsize-1)&&( (heap.asize & 0x80000000) != 0)  ) { // going to down_alloc
        heap.pop_back();
        rehash( hash_mag-1);
    } else{
    heap.pop_back();
    

    unsigned int i = hashpos( HashFunc::makeSeed(heap[ite].first.k) );
    unsigned int *j = hash + i;
    i = hash[i];
    for(;i != heap.size();i = heap[i].second){
        j = &(heap[i].second);
    }
    *j = ite;

    

    }
    }
  //  ExOp::show(heap);
}

LFHTEMP void myHashmap<Key,Data,HashFunc>::erase(const Key &k){
    
    unsigned int i = hashpos( HashFunc::makeSeed(k) );
    unsigned int *j = hash + i;
    i = hash[i];

    
    for(;i != 0xFFFFFFFF;i = heap[i].second){
        if (k == heap[i].first.k) break;
        j = &(heap[i].second);
        }
    
    if (i == 0xFFFFFFFF) return;
    
     
    *j = heap[i].second;
    swap_and_pop(i);

    } 
LFHTEMP void myHashmap<Key,Data,HashFunc>::erase_from_iterator(unsigned int ite){
    if (ite == 0xFFFFFFFF) return;
    unsigned int i = hashpos( HashFunc::makeSeed(heap[ite].first.k) );
    unsigned int *j = hash + i;

    while(*j != ite) j = &(heap[*j].second);
    *j = heap[*j].second;
    
    
    // swap and pop
    swap_and_pop(ite);
    }
LFHTEMP unsigned int myHashmap<Key,Data,HashFunc>::hashpos(unsigned int seed) const{return (seed & ((0xFFFFFFFF << hash_mag)^ 0xFFFFFFFF));}

LFHTEMP void myHashmap<Key,Data,HashFunc>::rehash(unsigned char _new_mag){
    if (hash_mag) delete[](hash);
    hash_mag = _new_mag;
 //   printf("rehashing to size %i\n", 1<<_new_mag); fflush(stdout);
    if (_new_mag == 0) return;
    hash = new unsigned int[1 << hash_mag];
    memset(hash, '\xFF',sizeof(unsigned int) << hash_mag);

    unsigned int i;
   for(i=0;i< heap.size();i++){
        unsigned int j = hashpos( HashFunc::makeSeed(heap[i].first.k) );
        heap[i].second = hash[j];
        hash[j] = i;
        
}
	
}

LFHTEMP void myHashmap<Key,Data,HashFunc>::show(FILE*f, int level)const{
    if (hash_mag == 0) { fprintf(f, "Un-initialized HashMap, nbitems=%i\n", heap.size());}
    else{
        fprintf(f, "HashMap hashsize=%i, nbitems=%i\n", 1 << hash_mag, heap.size());
        unsigned int i;
        unsigned int ttt =0;
        for(i=0;(i>> hash_mag) == 0;i++){
               printf("%i:",i);
               for(unsigned int j= hash[i]; j != 0xFFFFFFFF;j = heap[j].second) {printf("\t");ExOp::show(heap[j].first,f,2);ttt++;}
               printf("\n");
            }
        if (ttt != heap.size()) {printf("ill state! %i %i\n", heap.size(), ttt); exit(1);}
    }
}

	
#undef LFHTEMP
#define LFHTEMP template<class Key, class BackKey, class Data, class HashFunc, class BackHashFunc >
LFHTEMP void dualHashmap<Key,BackKey,Data,HashFunc,BackHashFunc>::set(const Key & key, const BackKey & Bkey){
    unsigned int j = find(key); 
    if (j == 0xFFFFFFFF){
        
        pair< KeyElem<pair <Key, BackKey> , Data> , pair<unsigned int,unsigned int> > new_buck;
        new_buck.first.k.first = key;
        new_buck.first.k.second = Bkey;
       	int rsize = (heap.asize & 0x7FFFFFFF);
        j = heap.size();
        if ((((heap.asize + 1)^(heap.asize))> rsize)&&((heap.asize & 0x80000000) == 0)) { // going to up_alloc 
            heap.push_back(new_buck);
            rehash(hash_mag + 1);
        }else{
            heap.push_back(new_buck);
            unsigned int i = hashpos( HashFunc::makeSeed(key) );
            heap[j].second.first = hash[i];
            hash[i] = j;

            i = hashpos( BackHashFunc::makeSeed(Bkey) ) | (1<<hash_mag);
            heap[j].second.second = hash[i];
            hash[i] = j;
        }

    }

}

LFHTEMP unsigned int dualHashmap<Key,BackKey,Data,HashFunc,BackHashFunc>::find(const Key & k )const{
    if (hash_mag == 0) return(0xFFFFFFFF);
    unsigned int j = hash[hashpos( HashFunc::makeSeed(k) )];
    for(;j != 0xFFFFFFFF;j = heap[j].second.first){
        if (k == heap[j].first.k.first) break;
    }
    return(j);
}

LFHTEMP unsigned int dualHashmap<Key,BackKey,Data,HashFunc,BackHashFunc>::back_find(const BackKey & b )const{
    if (hash_mag == 0) return(0xFFFFFFFF);
    unsigned int j = hash[hashpos( BackHashFunc::makeSeed(b) ) | (1 << hash_mag) ];
    for(;j != 0xFFFFFFFF;j = heap[j].second.second){
        if (b == heap[j].first.k.second) break;
    }
    return(j);
}

LFHTEMP void dualHashmap<Key,BackKey,Data,HashFunc,BackHashFunc>::swap_and_pop(unsigned int ite){
    // swap and pop
   // printf("pop swap %i, ,%i\n", ite, heap.size());fflush(stdout);
    int rsize = (heap.asize & 0x7FFFFFFF);
    if (ite == heap.size()-1){ // popped the last element of the array (no swap, check if array gets empty)
    if ((((heap.asize)^(heap.asize-1))> rsize-1)&&(( (heap.asize & 0x80000000) != 0)||(rsize == 1) )  ) {
        heap.pop_back();rehash((rsize == 1) ? 0 : hash_mag-1);
    }else{
        heap.pop_back();
    }

    }else{

    heap[ite] = heap[heap.size()-1];
    
    
    
    
    
    if ((((heap.asize)^(heap.asize-1))> rsize-1)&&( (heap.asize & 0x80000000) != 0)  ) { // going to down_alloc
        heap.pop_back();
        rehash( hash_mag-1);
    } else{
    heap.pop_back();
    

    unsigned int i = hashpos( HashFunc::makeSeed(heap[ite].first.k.first) );
    unsigned int *j = hash + i;
    i = hash[i];
    for(;i != heap.size();i = heap[i].second.first){
        j = &(heap[i].second.first);
    }
    *j = ite;

    i = hashpos( BackHashFunc::makeSeed(heap[ite].first.k.second) ) | (1 << hash_mag);
    j = hash + i;
    i = hash[i];
    for(;i != heap.size();i = heap[i].second.second){
        j = &(heap[i].second.second);
    }
    *j = ite;

    }
    }
  //  ExOp::show(heap);
}

LFHTEMP void dualHashmap<Key,BackKey,Data,HashFunc,BackHashFunc>::erase_from_iterator(unsigned int ite){
    if (ite == 0xFFFFFFFF) return;

    unsigned int *j = hash + hashpos( HashFunc::makeSeed(heap[ite].first.k.first) );

    
    while(*j != ite) j = &(heap[*j].second.first);
    *j = heap[*j].second.first;
    

    j = hash + (hashpos( BackHashFunc::makeSeed(heap[ite].first.k.second) ) | (1 << hash_mag));

    while(*j != ite) j = &(heap[*j].second.second);
    *j = heap[*j].second.second;

    // swap and pop
    swap_and_pop(ite);
    }

LFHTEMP void dualHashmap<Key,BackKey,Data,HashFunc,BackHashFunc>::rehash(unsigned char _new_mag){
    if (hash_mag) delete[](hash);
    hash_mag = _new_mag;
    printf("rehashing to size %i\n", 1<<_new_mag); fflush(stdout);
    if (_new_mag == 0) return;
    hash = new unsigned int[2 << hash_mag];
    memset(hash, '\xFF',sizeof(unsigned int) << (hash_mag +1));

    unsigned int i;
   for(i=0;i< heap.size();i++){
        unsigned int j = hashpos( HashFunc::makeSeed(heap[i].first.k.first) );
        heap[i].second.first = hash[j];
        hash[j] = i;
   
        j = hashpos(BackHashFunc::makeSeed(heap[i].first.k.second) ) | (1<< hash_mag);
        heap[i].second.second = hash[j];
        hash[j] = i;     
}
	
}

LFHTEMP void dualHashmap<Key,BackKey,Data,HashFunc,BackHashFunc>::show(FILE*f, int level)const{
    if (hash_mag == 0) { fprintf(f, "Un-initialized HashMap, nbitems=%i\n", heap.size());}
    else{
        fprintf(f, "HashMap hashsize=%i, nbitems=%i\n", 1 << hash_mag, heap.size());
        unsigned int i;
        unsigned int ttt =0;
        for(i=0;(i>> hash_mag) == 0;i++){
               printf("%i:",i);
               for(unsigned int j= hash[i]; j != 0xFFFFFFFF;j = heap[j].second.first) {printf("\t");ExOp::show(heap[j].first,f,2);ttt++;}
               printf("\n");
            }
        fprintf(f, "BackMap:\n");
        for(i=0;(i>> hash_mag) == 0;i++){
               printf("%i:",i);
               for(unsigned int j= hash[i | (1 << hash_mag) ]; j != 0xFFFFFFFF;j = heap[j].second.second) {printf("\t");ExOp::show(heap[j].first,f,2);ttt++;}
               printf("\n");
            }
        if (ttt != heap.size()) {printf("ill state! %i %i\n", heap.size(), ttt); exit(1);}
    }
}

LFHTEMP void dualHashmap<Key,BackKey,Data,HashFunc,BackHashFunc>::test(FILE*f){

Vector< KeyElem<Key,BackKey> > list;

KeyElem<Key,BackKey > tmpelem;
fprintf(f,"Dual hash map test, randomly insert and delete element, with objective size alternating between 8 and  2000 every 4096 iterations for million operations\n");
unsigned int i;
unsigned int r,j;
for(i=0;i<1000000;i++){
    ExOp::toRand(r); r >>= (i & 4096) ? 20: 28;
    if (r > list.size()){
        j=0;
        do{
        ExOp::toRand(tmpelem);
        r = find(tmpelem.k);
        if (r == 0xFFFFFFFF) r = back_find(tmpelem.d);
        if (j++ > 1000) {fprintf(f,"Could not find a used key-data pair to insert, exiting\n"); return;}
        }while (r != 0xFFFFFFFF);
        set(tmpelem.k, tmpelem.d); list.push_back(tmpelem);
        fprintf(f,"Insert :"); tmpelem.show(f,1);fprintf(f,"\n"); 
        
    }else{
        ExOp::toRand(r); j = r % list.size();
        r = find(list[j].k);
        fprintf(f,"Delete :"); list[j].show(f,1);fprintf(f," at %i\n", r);
        if (r == 0xFFFFFFFF) exit(1);
        if (back_deref(r) != list[j].d) exit(1);
        if (front_deref(r) != list[j].k) exit(1);
        if (r != back_find(list[j].d)) exit(1);
        erase_from_iterator(r);
        list.pop_swap(j);
        
    }

}

}


 // namespace end
