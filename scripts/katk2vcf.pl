#!/usr/bin/perl

my $chr_dir = "/storage9/db/human_GRCh38/data/chr/";

@chrs = ('MT','X','Y','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22');
$ref = "GRCh38";
$start = time();
foreach $chr (@chrs) {
   $chr_seq = "";
   $chr_file = "${chr_dir}${chr}.fa";
   open (CHR,"$chr_file");
   while(<CHR>){
      chomp;
      next if /^\s*$/;
      next if /^\s*>/;
      $chr_seq .= $_;
   }
   $hg38{$chr} = $chr_seq;
   close(CHR);
}
                                             

$rida = 0;
open CALL,"$ARGV[0]" or die;
@jarjestus = @muutused = @muutused_nc = @posit = ();
%call = %vcf = %sim = %nc = ();
$pikk = 0;

$ref_tulp = 3;
$call_tulp = 5;
$tyyp_tulp = 6;
$eelmine_tulp = 9;
$jrk = -1;

while(<CALL>){
  chomp;
  $rida++;
  @tmp0 = split(/\t/);
  $posit[$rida] = $tmp0[1];
  $jarjestus[$rida] = $tmp0[$ref_tulp];
  $key0 = $tmp0[0]; $key0 .= ":"; $key0 .= $tmp0[1];
  $key1 = $tmp0[0]; $key1 .= ":"; $key1 .= $tmp0[1]-1;
  if($tmp0[5] eq "NC"){ $nc{$key0} = "NC"; push @muutused_nc, $key0; next;}
                 
  if($pikk == 1 && $posit[$rida] - $posit[$rida-1] > 1){
      @lahti = split(/\:/,$voti);     
      $lahti[1] = $lahti[1] -1;
      if($taht eq "I"){
         $ref_I_2 = $true_I_2 = $ajut_2 = $ajut_mut_2 = $ajut_pos = "";
         $mut_pikk_2 = length($mutat_2);
         $ref_I_2 = substr($hg38{$lahti[0]},$lahti[1]-50,100) if length($mutat_2) > 1;
         $true_I_2 = substr($hg38{$lahti[0]},$lahti[1]-50,51) if length($mutat_2) > 1;
         $true_I_2 .= substr($mutat_2,1) if length($mutat_2) > 1;
         $true_I_2 .= substr($hg38{$lahti[0]},$lahti[1]+1,49) if length($mutat_2) > 1;
         for($r = 0; $r < 50; $r++){
            $ajut_2 = substr($true_I_2,0,50-$r);
            $ajut_2 .= substr($true_I_2,-50-$r);
            if($ajut_2 eq $ref_I_2 && length($mutat_2) > 1) { 
               $ajut_pos = $lahti[1]-$r;
               $mutat_1 = substr($true_I_2,49-$r,1) if $het == 1; 
               $mutat_2 = substr($true_I_2,49-$r,$mut_pikk_2);
               $mutat_1 = $mutat_2 if $het == 0; 
               $voti = $tmp0[0]; 
               $voti .= ":"; 
               $voti .= $ajut_pos;
            }
         }
      }
      
      if($taht eq "D"){
         $ref_D_1 = $ajut_1 = $ajut_pos = "";
         $ref_D_1 = substr($hg38{$lahti[0]},$lahti[1]-50,50) if length($mutat_1) > 1;
         $ref_D_1 .= substr($hg38{$lahti[0]},$lahti[1]+length($mutat_1)-1,50) if length($mutat_1) > 1;   
         for($r = 0; $r < 50; $r++){
            $ajut_1 = substr($hg38{$lahti[0]},$lahti[1]-50,50-$r) if length($mutat_1) > 1;
            $ajut_1 .= substr($hg38{$lahti[0]},$lahti[1]+length($mutat_1)-1-$r,50+$r) if length($mutat_1) > 1;
            if($ajut_1 eq $ref_D_1) { 
               $ajut_pos = $lahti[1]-$r; 
               $mutat_2 = lc(substr($hg38{$lahti[0]},$lahti[1]-$r-1,1)) if $het == 1;
               $mutat_1 = lc(substr($hg38{$lahti[0]},$lahti[1]-$r-1,length($mutat_1)));
               $mutat_2 = $mutat_1 if $het == 0; 
               $voti = $tmp0[0]; 
               $voti .= ":"; 
               $voti .= $ajut_pos;
            }
         }
         $mutat_1 = ucfirst($mutat_1);
         $mutat_2 = ucfirst($mutat_2);                                                                                                                    
      }                                                


      $mutat = $mutat_1; $mutat .= "/"; $mutat .= $mutat_2;
      $jrk++;
      $callid[$jrk] = $voti;
      $call{$voti} = $mutat;
      $tyyp{$voti} = $taht;
      $muut{$voti} = 1 unless exists $muut{$voti};
      $pikk = 0;
      $taht = "";
      $het = 0;
  }

  # insertion "I"
   
  if ($tmp0[$tyyp_tulp] eq "I" && $pikk == 0){
      $voti = $key0;
      $mutat_1 = $mutat_2 = substr($hg38{$tmp0[0]},$tmp0[1]-1,1);
      if(substr($tmp0[$call_tulp],0,1) ne substr($tmp0[$call_tulp],1,1)){ 
         $mutat_2 .= substr($tmp0[$call_tulp],1,1) if substr($tmp0[$call_tulp],0,1) eq "-";
         $mutat_2 .= substr($tmp0[$call_tulp],0,1) if substr($tmp0[$call_tulp],1,1) eq "-";
         $het = 1;
      }
      if(substr($tmp0[$call_tulp],0,1) eq substr($tmp0[$call_tulp],1,1)){
         $mutat_1 .= substr($tmp0[$call_tulp],0,1);
         $mutat_2 .= substr($tmp0[$call_tulp],1,1);
      }
      $pikk = 1;
      $taht = "I";
      next;
  }
  
  if ($tmp0[$tyyp_tulp] eq "I" && $pikk == 1){
      if(substr($tmp0[$call_tulp],0,1) ne substr($tmp0[$call_tulp],1,1)){
         $mutat_2 .= substr($tmp0[$call_tulp],1,1) if substr($tmp0[$call_tulp],0,1) eq "-";
         $mutat_2 .= substr($tmp0[$call_tulp],0,1) if substr($tmp0[$call_tulp],1,1) eq "-";
      }
      if(substr($tmp0[$call_tulp],0,1) eq substr($tmp0[$call_tulp],1,1)){
         $mutat_1 .= substr($tmp0[$call_tulp],0,1);
         $mutat_2 .= substr($tmp0[$call_tulp],0,1);
      }
      next;
  }
  
  # insertion "D"  
  
  if ($tmp0[$tyyp_tulp] eq "D" && $pikk == 0){
     $voti = $key0;
     $mutat_1 = $mutat_2 = substr($hg38{$tmp0[0]},$tmp0[1]-1,1);
     if(substr($tmp0[$call_tulp],0,1) ne substr($tmp0[$call_tulp],1,1)){
        $mutat_2 .= lc(substr($tmp0[$call_tulp],1,1)) if substr($tmp0[$call_tulp],0,1) eq "-";
        $mutat_1 .= lc(substr($tmp0[$call_tulp],0,1)) if substr($tmp0[$call_tulp],1,1) eq "-";
        $het = 1;
     }
     if(substr($tmp0[$call_tulp],0,1) eq substr($tmp0[$call_tulp],1,1)){
        $mutat_1 .= lc(substr($tmp0[$call_tulp],0,1));
        $mutat_2 .= lc(substr($tmp0[$call_tulp],0,1));
     }
     $pikk = 1;
     $taht = "D";
     next;
  }
  
  if ($tmp0[$tyyp_tulp] eq "D" && $pikk == 1){
     if(substr($tmp0[$call_tulp],0,1) ne substr($tmp0[$call_tulp],1,1)){
        $mutat_2 .= lc(substr($tmp0[$call_tulp],1,1)) if substr($tmp0[$call_tulp],0,1) eq "-";
        $mutat_1 .= lc(substr($tmp0[$call_tulp],0,1)) if substr($tmp0[$call_tulp],1,1) eq "-";
     }
     if(substr($tmp0[$call_tulp],0,1) eq substr($tmp0[$call_tulp],1,1)){
        $mutat_1 .= lc(substr($tmp0[$call_tulp],0,1));
        $mutat_2 .= lc(substr($tmp0[$call_tulp],0,1));
     }
     next;
  }    
  
  # substitution "S"
  
  if ($tmp0[$tyyp_tulp] eq "S"){
     if(substr($tmp0[$call_tulp],0,1) ne substr($tmp0[$call_tulp],1,1)){ 
        $mutat_1 = $tmp0[$ref_tulp];
        $mutat_2 = substr($tmp0[$call_tulp],0,1) if substr($tmp0[$call_tulp],1,1) eq $tmp0[$ref_tulp];
        $mutat_2 = substr($tmp0[$call_tulp],1,1) if substr($tmp0[$call_tulp],0,1) eq $tmp0[$ref_tulp];
     }
     if(substr($tmp0[$call_tulp],0,1) eq substr($tmp0[$call_tulp],1,1)){ $mutat_1 = substr($tmp0[$call_tulp],0,1); $mutat_2 = substr($tmp0[$call_tulp],1,1);}
     $mutat = $mutat_1; $mutat .= "/"; $mutat .= $mutat_2;
     $jrk++;
     $callid[$jrk] = $key0;
     $call{$key0} = $mutat;
     $tyyp{$key0} = $tmp0[$tyyp_tulp];
     $muut{$key0} = 1 unless exists $muut{$key0};
     next;
  }          
}
close CALL;
                                                                              
print "##fileformat=VCFv4.0\n";
print "##fileDate=\n";
print "##source=KATKtools\n";
print "##reference=$ref\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT$sample\n";

for($l = 0; $l < $jrk; $l++){
   $muutus = $callid[$l];
   @asukoht = split(/\:/,$muutus);
   $call{$muutus} = uc($call{$muutus});
   @genotyybid = split(/\//,$call{$muutus});
   $nuc = substr($hg38{$asukoht[0]},$asukoht[1]-1,1);
   $call{$muutus} = $nc{$muutus} if exists $nc{$muutus};
   if ($tyyp{$muutus} eq "I" || $tyyp{$muutus} eq "S"){
      print"$asukoht[0]\t$asukoht[1]\t.\t$nuc\t$genotyybid[1]\t.\tPASS\t$tyyp{$muutus}\tGT\t";
      print"0" if $genotyybid[0] eq $nuc;
      print"1" if $genotyybid[0] eq $genotyybid[1];
      print"/";
      print"1\n";
      next;
   }
   if ($tyyp{$muutus} eq "D"){
      print"$asukoht[0]\t$asukoht[1]\t.\t$genotyybid[0]\t$nuc\t.\tPASS\t$tyyp{$muutus}\tGT\t"; 
      print"0";
      print"/";
      print"0\n" if $genotyybid[0] eq $genotyybid[1];
      print"1\n" if $genotyybid[0] ne $genotyybid[1];
   }
}
