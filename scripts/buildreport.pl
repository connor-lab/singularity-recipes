#!/usr/bin/env perl

use RTF::Writer;
use Getopt::Long;
use JSON;
use Data::Dumper;
use File::Slurp;
use feature qw/ say /;
use DateTime;

my $dt = DateTime->today(time_zone=>'local');
my $date = $dt->dmy('-');

#my $format = DateTime::Format::Strptime->new( pattern => '%d%m%Y' );
#$date = $format->format_datetime($date);
#use Data::Diver qw( Dive );

# so we need to, query sierra, get the reply and parse that json and ouput a PDF

Getopt::Long::Configure ('bundling');
GetOptions ('i=s' => \$sequence,
            'n=s' => \$sample,
            'o=s' => \$operator,
            'l=s' => \$lab,
            'h' => \$help) or help_msg();

$log_file = $sample.'.log';
$jsonfile = $sample.".json";
my $cmd = "sierrapy fasta ".$sequence." > ".$jsonfile;
# log execution
open my $logfh, '>>', $log_file;
print $logfh "\n### ", $cmd, "\n\n";
close $logfh;
# run it
$cmd .= " 2>> $log_file";
print "Running: $cmd\n";
system($cmd)==0 or err("Error running command, check $outdir/$log_file");

my $file = read_file($jsonfile);

my $perljson  = decode_json $file;


my $obj = from_json( $file );

print "from JSON\n";

@prtext;
@rttext;
@nrttext;
@intext;

%rtvals;
%prvals;
%invals;
%nrvals;

$rtrange ="";
$prrange ="";
$inrange ="";

$hivdbv ="";
$hivdbd ="";
# add mutation hash

# add comment hash for each

# 0 = major ; 1 = minor ; 2 = other ; 3 = comments ; 4 = resistances

for($i = 0 ; $i <5 ; $i++)
{
$prtext[$i] =" ";
$rttext[$i] =" ";
$nrttext[$i] =" ";
$intext[$i] =" ";
}


my $gene = "";
my $drugres = "";
my $drugscores = "";
for my $vals (@$perljson) {
    $subtypetxt = $subtypetxt ." ".$vals->{'subtypeText'};
    $drugres = $vals->{'drugResistance'};
    for my $drugresvs (@$drugres) {
      $hivdbv = $drugresvs->{'version'}->{'text'};
      $hivdbd = $drugresvs->{'version'}->{'publishDate'};

      $gene = $drugresvs->{'gene'}->{'name'};
      $drugscores = $drugresvs->{'drugScores'};
      for my $inDrugScores (@$drugscores)
      {
	print $gene "\n";
	if($gene =~ m/RT/)
	{
		if($inDrugScores->{'drugClass'}->{'name'} =~ m/NNRTI/)
		{
      $nrttext[4] = $nrttext[4].$inDrugScores->{'drug'}->{'name'}."\t".$inDrugScores->{'text'}."\n";
      #print $inDrugScores->{'text'}."\n";
      if($inDrugScores->{'text'} =~ m/Susceptible/)
        {

           }
      elsif(($inDrugScores->{'text'} =~ m/Low-Level Resistance/)||($inDrugScores->{'text'} =~ m/Intermediate Resistance/)||($inDrugScores->{'text'} =~ m/High-Level Resistance/))
        {
            # 0 = major ; 1 = minor ; 2 = other ; 3 = comments ; 4 = resistances
            $partialScores = $inDrugScores->{'partialScores'};
            for my $partials (@$partialScores) {
            $keyval = $partials->{'mutations'}->[0]->{'text'};
            $keycontent = $partials->{'mutations'}->[0]->{'comments'}->[0]->{'text'};
            #if (!exists($nrvals{$keyval})) {
            if (index($nrttext[0], $keyval) == -1) {
              $nrttext[0] = $nrttext[0].$partials->{'mutations'}->[0]->{'text'}.", ";
              $nrttext[3] = $nrttext[3].$partials->{'mutations'}->[0]->{'comments'}->[0]->{'text'}."\n";
              }
          }
        }
      }
    elsif($inDrugScores->{'drugClass'}->{'name'} =~ m/NRTI/)
		{
      $rttext[4] = $rttext[4].$inDrugScores->{'drug'}->{'name'}."\t".$inDrugScores->{'text'}."\n";
      #print $inDrugScores->{'text'}."\n";
      if($inDrugScores->{'text'} =~ m/Susceptible/)
        {
			       #$rttext[4] = $rttext[4].$inDrugScores->{'drug'}->{'name'}."\t".$inDrugScores->{'text'}."<br>";
        }
      elsif(($inDrugScores->{'text'} =~ m/Low-Level Resistance/)||($inDrugScores->{'text'} =~ m/Intermediate Resistance/)||($inDrugScores->{'text'} =~ m/High-Level Resistance/))
        {
          # 0 = major ; 1 = minor ; 2 = other ; 3 = comments ; 4 = resistances
          $partialScores = $inDrugScores->{'partialScores'};
          for my $partials (@$partialScores) {
            $keyval = $partials->{'mutations'}->[0]->{'text'};
            $keycontent = $partials->{'mutations'}->[0]->{'comments'}->[0]->{'text'};
            if (index($rttext[0], $keyval) == -1) {
              $rttext[0] = $rttext[0].$partials->{'mutations'}->[0]->{'text'}.", ";
              $rttext[3] = $rttext[3].$partials->{'mutations'}->[0]->{'comments'}->[0]->{'text'}."\n";
            }

           }
        }
		}
	}
	elsif($gene =~ m/PR/)
	{
    $prtext[4] = $prtext[4].$inDrugScores->{'drug'}->{'name'}."\t".$inDrugScores->{'text'}."\n";
      if($inDrugScores->{'text'} =~ m/Susceptible/)
        {

             }
      elsif(($inDrugScores->{'text'} =~ m/Low-Level Resistance/)||($inDrugScores->{'text'} =~ m/Intermediate Resistance/)||($inDrugScores->{'text'} =~ m/High-Level Resistance/))
        {
          $partialScores = $inDrugScores->{'partialScores'};
          for my $partials (@$partialScores) {
          $keyval = $partials->{'mutations'}->[0]->{'text'};
          $keycontent = $partials->{'mutations'}->[0]->{'comments'}->[0]->{'text'};
          if($partials->{'mutations'}->[0]->{'primaryType'} =~ m/Major/)
          {
            if (index($prtext[0], $keyval) == -1) {
              $prtext[0] = $prtext[0].$partials->{'mutations'}->[0]->{'text'}.", ";
              $prtext[3] = $prtext[3].$partials->{'mutations'}->[0]->{'comments'}->[0]->{'text'}."\n";
              }
            #$prtext[0] = $prtext[0].$inDrugScores->{'partialScores'}->[0]->{'mutations'}->[0]->{'text'}.", ";
          }
          if($partials->{'mutations'}->[0]->{'primaryType'} =~ m/Minor/)
          {
            if (index($prtext[1], $keyval) == -1) {
              $prtext[1] = $prtext[1].$partials->{'mutations'}->[0]->{'text'}.", ";
              $prtext[3] = $prtext[3].$partials->{'mutations'}->[0]->{'comments'}->[0]->{'text'}."\n";
              }
            #$prtext[1] = $prtext[1].$partials->{'mutations'}->[0]->{'text'}.", ";
          }
            #$prtext[3] = $prtext[3].$partials->{'mutations'}->[0]->{'comments'}->[0]->{'text'}.", ";
          }
        }

    }
	elsif($gene =~ m/IN/)
		{
      $intext[4] = $intext[4].$inDrugScores->{'drug'}->{'name'}."\t".$inDrugScores->{'text'}."\n";
      if($inDrugScores->{'text'} =~ m/Susceptible/)
      {

			}
      elsif(($inDrugScores->{'text'} =~ m/Low-Level Resistance/)||($inDrugScores->{'text'} =~ m/Intermediate Resistance/)||($inDrugScores->{'text'} =~ m/High-Level Resistance/))
      {
        $partialScores = $inDrugScores->{'partialScores'};
        for my $partials (@$partialScores) {
        $keyval = $partials->{'mutations'}->[0]->{'text'};
        $keycontent = $partials->{'mutations'}->[0]->{'comments'}->[0]->{'text'};

        if($partials->{'mutations'}->[0]->{'primaryType'} =~ m/Major/)
        {
          if (index($intext[0], $keyval) == -1) {
            $intext[0] = $intext[0].$partials->{'mutations'}->[0]->{'text'}.", ";
            $intext[3] = $intext[3].$partials->{'mutations'}->[0]->{'comments'}->[0]->{'text'}."\n";
          }
          #$intext[0] = $intext[0].$partials->{'mutations'}->[0]->{'text'}.", ";
        }
        if($partials->{'mutations'}->[0]->{'primaryType'} =~ m/Minor/)
        {
          if (index($intext[1], $keyval) == -1) {
            $intext[1] = $intext[1].$partials->{'mutations'}->[0]->{'text'}.", ";
            $intext[3] = $intext[3].$partials->{'mutations'}->[0]->{'comments'}->[0]->{'text'}."\n";
          }
          #$intext[1] = $intext[1].$partials->[0]->{'mutations'}->[0]->{'text'}.", ";
        }
          #$intext[3] = $intext[3].$partials->[0]->{'mutations'}->[0]->{'comments'}->[0]->{'text'}.", ";
        }
      }

		}

        #print "Drug: ".$inDrugScores->{'drug'}->{'name'}."\n Class: ".$inDrugScores->{'drugClass'}->{'name'}."\n Score: ".$inDrugScores->{'score'}."\n Status: ".$inDrugScores->{'text'}."\n";

      }

    }
    $alignedGeneSequences = $vals->{'alignedGeneSequences'};
    for my $geneseqs(@$alignedGeneSequences)
    {
      $gene = $geneseqs->{'gene'}->{'name'};
      if($gene =~ m/PR/)
      {
          $prrange =  $geneseqs->{'firstAA'}."-".$geneseqs->{'lastAA'};
          $mutations = $geneseqs->{'mutations'};
          for my $partials (@$mutations) {
          $keyval = $partials->{'consensus'}.$partials->{'position'}.$partials->{'AAs'};
          #print $keyval."\n";
          if ((index($prtext[0], $keyval) == -1) || (index($prtext[1], $keyval) == -1)) {
            $prtext[2] = $prtext[2].$keyval.", ";
          }
        }
      }
      elsif($gene =~ m/RT/)
      {
          $rtrange =  $geneseqs->{'firstAA'}."-".$geneseqs->{'lastAA'};
          $mutations = $geneseqs->{'mutations'};
          for my $partials (@$mutations) {
          $keyval = $partials->{'consensus'}.$partials->{'position'}.$partials->{'AAs'};
          #print $keyval."\n";
          if ((index($rttext[0], $keyval) == -1) || (index($nrttext[0], $keyval) == -1)) {
            $rttext[2] = $rttext[2].$keyval.", ";
          }
        }
      }
      elsif($gene =~ m/IN/)
      {
          $inrange =  $geneseqs->{'firstAA'}."-".$geneseqs->{'lastAA'};
          $mutations = $geneseqs->{'mutations'};
          for my $partials (@$mutations) {
          $keyval = $partials->{'consensus'}.$partials->{'position'}.$partials->{'AAs'};
          #print $keyval."\n";
          if ((index($intext[0], $keyval) == -1) || (index($intext[1], $keyval) == -1)) {
            $intext[2] = $intext[2].$keyval.", ";
          }
        }

      }
    }
}

$subtypetxt =~ s/NA//g;
$subtypetxt =~ s/ //g;
my $rtf = RTF::Writer->new_to_file($sample.".rtf");
$rtf->prolog('title' => 'PHW HIV Report for sample'.$sample, 'fonts' => ["Verdana"], 'colors' => [undef, [188,188,188], [151,151,151]]);

#$rtf->number_pages;
$rtf->paragraph(
\'\qc\fs40\b\i',
$rtf->image( 'filename' => "/usr/local/bin/phw.jpg", ),
"\nHIV Resistance Genotyping Report\n"
);
$rtf->paragraph(
\'\qc\fs28\b\i',
"PHW Pathogen Genomics Unit (PenGU), Public Health Wales Microbiology, University Hospital Wales, Cardiff\n"
);

my $tablehead = RTF::Writer::TableRowDecl->new('widths' => [9000]);
$rtf->row($tablehead, [\'\b',"Sender\'s details"]);
my $table1 = RTF::Writer::TableRowDecl->new('widths' => [2250,2250,2250,2250]);
$rtf->row($table1, "Requesting Clinician","", "Sender\'s Lab Number", "");
$rtf->row($table1, "Location", "", "Sample Date", "");
my $tablehead2 = RTF::Writer::TableRowDecl->new('widths' => [9000]);
$rtf->row($tablehead2, [\'\b',"Patient and Sample Information"]);
my $table2 = RTF::Writer::TableRowDecl->new('widths' => [2250,2250,2250,2250]);
$rtf->row($table2, "Name", "", "PHW Episode Number", $sample);
$rtf->row($table2, "Patient ID number", "", "Date Sample Received","");
$rtf->row($table2, "DOB","","Sample Type","");
$rtf->row($table2, "Date of Report",$date,"Current Treatment","");
$rtf->row($table2, "Sample Viral Load","","Previous Treatment","");

$rtf->paragraph(\'\fs30',"\n");

my $prtablehead = RTF::Writer::TableRowDecl->new('widths' => [9000]);
$rtf->row($prtablehead, [\'\b','Protease (PR) codons analysed: '.$prrange]);
my $prtable = RTF::Writer::TableRowDecl->new('widths' => [2000,7000]);
$rtf->row($prtable, [\'\b','PI Major Resistance Mutations'],$prtext[0]);
$rtf->row($prtable, [\'\b','PI Minor Resistance Mutations'],$prtext[1]);
$rtf->row($prtable, [\'\b','Other Mutations'],$prtext[2]);
$rtf->row($prtable, [\'\b','Protease Inhibitors'], $prtext[4]);

$rtf->paragraph(\'\fs30',"\n");

my $rttablehead = RTF::Writer::TableRowDecl->new('widths' => [9000]);
$rtf->row($rttablehead, [\'\b','Reverse Transcriptase codons analysed: '.$rtrange]);
my $rttable = RTF::Writer::TableRowDecl->new('widths' => [4500,4500]);
$rtf->row($rttable, [\'\b','NRTI Resistance Mutations'],$rttext[0]);
$rtf->row($rttable, [\'\b','NNRTI Resistance Mutations'],$nrttext[0]);
$rtf->row($rttable, [\'\b','Other Mutations'],$rttext[2]);
$rtf->row($rttable, [\'\b','Nuceloside RTI'],[\'\b\cb1','Non-nucleoside RTI']);
$rtf->row($rttable, $rttext[4], $nrttext[4]);

$rtf->paragraph(\'\fs30',"\n");

my $inttablehead = RTF::Writer::TableRowDecl->new('widths' => [9000]);
$rtf->row($inttablehead, [\'\b','Integrase codons analysed: '.$inrange]);
my $inttable= RTF::Writer::TableRowDecl->new('widths' => [2000,7000]);
$rtf->row($inttable, [\'\b','IN Major Resistance Mutations'],$intext[0]);
$rtf->row($inttable, [\'\b','IN Minor Resistance Mutations'],$intext[1]);
$rtf->row($inttable, [\'\b','Other Mutations'],$intext[2]);
$rtf->row($inttable, [\'\b','Integrase strand transfer inhibitors'],$intext[4]);

$rtf->paragraph(\'\fs30',"\n");

my $gendets=RTF::Writer::TableRowDecl->new('widths' => [9000]);
$rtf->row($gendets, [\'\b','HIV Subtype: '.$subtypetxt]);
$rtf->row($gendets, [\'\b','PR Comments']);
$rtf->row($gendets, $prtext[3]);
$rtf->row($gendets, [\'\b','RT Comments']);
my $gendets2=RTF::Writer::TableRowDecl->new('widths' => [1100,3400,1100,3400]);
$rtf->row($gendets2,'NRTI:',$rttext[3],'NNRTI:',$nrttext[3]);
my $gendets3=RTF::Writer::TableRowDecl->new('widths' => [9000]);
$rtf->row($gendets3, [\'\b','IN Comments']);
$rtf->row($gendets3, $intext[3]);

$rtf->paragraph(
\'\qc\fs20\i',           
'
The results presented here are based upon the HIVdb Drug resistance database version '.$hivdbv.' released on the '.$hivdbd.'. The results were generated using the sierrapy web client on '.$date."\n"
);

my $signtab=RTF::Writer::TableRowDecl->new('widths' => [1250,3250,1250,3250]);
$rtf->row($signtab, "Signed:\n","","Name:","");
$rtf->row($signtab, "Position:\n","","Date:","");

$rtf->paragraph(\'\qc\fs20\i\b','
For technical enquiries please contact Dr Sally Corden (02920745226) sally.corden@wales.nhs.uk or Joanne Watkins (02920742046) joanne.watkins@wales.nhs.uk
For clinical enquires please contact Dr Nicola Price (02920742178) nicola.price@wales.nhs.uk or Dr Matthijs Backx (02920742166) matthijs.backx@wales.nhs.uk');


sub help_msg
  {
    print "buildreport.pl Options;
    -i input file containing HIV consensus sequence,
    -n sample name
    -o operator name,
    -l laboratory identifier,
    -h show this message\n";
  }
