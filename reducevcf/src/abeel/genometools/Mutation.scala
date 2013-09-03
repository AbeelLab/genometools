package abeel.genometools

/**
 * @author Thomas Abeel
 */
class Mutation(line: String) extends VCFLine(line) {

  def this(vcf:VCFLine)=this(vcf.line)
  
  lazy val identifier: String = {
    variation + "." + pos + "." + ref+"."+alt

  }

}