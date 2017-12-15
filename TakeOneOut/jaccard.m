function SimilarityIndex = jaccard(setA, setB)
%%%set A must be orginal


%%compute intersection 
inter = size(intersect(setA,setB),1);

%%union
uni = size(union(setA,setB),1);
%inter = setA == setB;
%den = size(union(setA,setB),1);
%num = sum(inter(:));
%den = sum(uni(:));

SimilarityIndex = inter/uni;
%SimilarityIndex = num/den;
%dis = pdist2(setA,setB,'jaccard');
%[a,b] = size(dis);
%val = sum(dis(:))/(a*b);
%SimilarityIndex = 1 - val;
end