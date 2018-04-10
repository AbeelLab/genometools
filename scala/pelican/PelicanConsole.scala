package pelican

object PelicanConsole {

  def main(args: Array[String]): Unit = {

    if (args.length == 0) {

      listInstructions
    } else {
      args(0) match {
        case "list" => listInstructions
        case "help" => listInstructions
        case "ambiguity" => AggregateAmbiguous.main(args.drop(1))
        case "picard" => PicardLocations.main(args.drop(1))
         case "pilon-metrics" => PilonMetrics.main(args.drop(1))
        case "snpmatrix" => PairwiseSNPdistance.main(args.drop(1))
        case "snpmatrix-distribution"=>PairwiseSNPdistribution.main(args.drop(1))
        case "suspicious" => FlagSuspectSamples.main(args.drop(1))
        case "variant2gene" => LSP2genes.main(args.drop(1))
        case "merge-projects" => MergeProjects.main(args.drop(1))
        case "extract-duplications"=>ExtractDuplications.main(args.drop(1))
         case "list-tree" => TreeOrder.main(args.drop(1))
        case _ => listInstructions
      }
    }

  }

  def listInstructions() {
    println("Usage:java -jar pelican.jar [instrucFtion] [instruction options...]")
    println("Instructions:")
    println("\tambiguity              Aggregate ambiguity information from VCF files")
    println("\tpicard                 Extract seq libraries locations from Picard")
    println("\tpilon-metrics          Program to extract metrics from all Pilon runs")
    println("\tsnpmatrix              Create pairwise SNP distance matrix")
    println("\tsnpmatrix-distribution Create distribution plot from SNP matrix")
    println("\tsuspicious             Identify samples that should be excluded based on their CPT, "+
        "ambiguity rate and large insertion rate.")
    println("\tvariant2gene           Map variants from VCF files to genes or other coordinate based features")
    println("\tmerge-projects         Merge several projects together.")
    println("\textract-duplications   Extract information about large duplications from Pilon runs")
    println("\tlist-tree              List all nodes in the tree in the order they will appear in the figure.")

    println("\tlist\t\tShow help")
    println("\thelp\t\tShow help")

  }

}
