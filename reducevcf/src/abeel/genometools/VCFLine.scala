package abeel.genometools

sealed trait SNP {

}
sealed trait Variation {
  override def toString():String={
   strType
  }
  lazy val strType = this.getClass().getSimpleName()

}

sealed trait LP;

class Match extends Variation with SNP 

class Substitution extends Variation
case class SingleSubstitution extends Substitution with SNP 
case class LongSubstitution extends Substitution with LP
//case class ComplexSubstitution extends Variation

class Deletion extends Variation
case class SingleDeletion extends Deletion with SNP 
case class LongDeletion extends Deletion with LP

class Insertion extends Variation
case class SingleInsertion extends Insertion with SNP
case class LongInsertion extends Insertion with LP

class VCFLine(val line: String) {

  private lazy val arr = line.split("\t")

  override def toString() = line //arr.mkString("\t")

  lazy val zeroPos = arr(1).toInt - 1

  lazy val pos = arr(1).toInt

  def ref = arr(3)
  def alt = arr(4)

  def refGenome = arr(0)

  lazy val refLength = ref.length()
  lazy val altLength = alt.length()

  lazy val variation: Variation = {
    if(alt.equals(".")){
      new Match
    }else if (ref.length() == alt.length()) {
      if (ref.equals(alt))
        new Match
      else {
        if (ref.length() == 1)
          new SingleSubstitution
        else
          new LongSubstitution
      }
    } else {
      assume(ref.length() > 0 && alt.length() > 0)
      if (ref.length() == 1 || alt.length() == 1) {
        val diff = ref.length() - alt.length()
        //if (ref.length() > alt.length())
        if (diff > 1)
          new LongDeletion
        else if (diff > 0)
          new SingleDeletion
        else if (diff < -1)
          new LongInsertion
        else if (diff < 0)
          new SingleInsertion
        else
          throw new RuntimeException("This is not supposed to happen!")
      } else {
        new LongSubstitution
      }

    }
  }

  lazy val score= if(arr(5)(0)=='.')0 else arr(5).toDouble
  
  lazy val pass = filter.equals("PASS")

  def filter = if (arr.length > 6) arr(6) else "."

}