//package ara

package abeel.genometools.vcf

/** Unapply methods for reading mutations from reduced VCFs*/

abstract class Mutation {
  type Mutation
  type SNP <: Mutation
  type Insertion <: Mutation
  type Deletion <: Mutation
  override def toString(): String 
  //def compare(that: Object): Int 
}

object Mutation {
  class Line(s: String) {
    def isSNP: Boolean = s match {
      case SNP(r, c, a,chr) => true
      case _ => false
    }
    
    def isValid: Boolean = s match {
      case SNP(r, c, a,chr) => true
      case Insertion(r, c, a,chr) => true
      case Deletion(r, c, a,chr) => true
      case MNP(r, c, a,chr) => true
      case _ => false
    }
    
    def invalidSite: Boolean = {
      val arr = s.split("\t")
      if (arr(4) == "." && arr(6) == "PASS") true
      else false
    }
    
  }
  implicit def seqtoBool(s: String) = new Line(s)

  case class SNP(val ref: String, val coordination: Int, val alt: String,val chr:String) extends Mutation {//} with Ordered[SNP] {
    override def toString(): String = ref + coordination + alt+chr
    //def compare(that: SNP): Int = this.coordination compare that.coordination
  }

  object SNP {
    def unapply(line: String): Option[(String, Int, String,String)] = {
      val arr = line.mkString.split("\t")
      val filter = arr(6)
      val ref = arr(3)
      val alt = arr(4)
      val coordination = arr(1).toInt
      val chromosome=arr(0)
      val nucleotides = Array[String]("A", "C", "T", "G")
      if (filter == "PASS" && nucleotides.contains(ref) && nucleotides.contains(alt)) Some((ref, coordination, alt,chromosome))
      else None
    }
  }

  case class Insertion(val ref: String, val coordination: Int, val alt: String,val chr:String) extends Mutation {//}with Ordered[Insertion] {
    override def toString(): String = ref + coordination + alt+chr
    //def compare(that: Insertion): Int = this.coordination compare that.coordination
  }

  object Insertion {
    def unapply(line: String): Option[(String, Int, String,String)] = {
      val arr = line.mkString.split("\t")
      val filter = arr(6)
      val ref = arr(3)
      val alt = arr(4)
      val coordination = arr(1).toInt
      val chromosome=arr(0)
      if (filter == "PASS" && ref.size < alt.size) Some((ref, coordination, alt,chromosome))
      else None
    }
  }

  case class Deletion(val ref: String, val coordination: Int, val alt: String,val chr:String) extends Mutation {//}with Ordered[Deletion] {
    override def toString(): String = ref + coordination + alt+chr
    //def compare(that: Deletion): Int = this.coordination compare that.coordination
  }

  object Deletion {
    def unapply(line: String): Option[(String, Int, String,String)] = {
      val arr = line.mkString.split("\t")
      val filter = arr(6)
      val ref = arr(3)
      val alt = arr(4)
      val coordination = arr(1).toInt
      val chromosome=arr(0)
      if (filter == "PASS" && ref.size > alt.size) Some((ref, coordination, alt,chromosome))
      else None
    }
  }
  
  case class MNP(val ref: String, val coordination: Int, val alt: String,val chr:String) extends Mutation {//}with Ordered[Deletion] {
    override def toString(): String = ref + coordination + alt+chr
    //def compare(that: Deletion): Int = this.coordination compare that.coordination
  }

  object MNP {
    def unapply(line: String): Option[(String, Int, String,String)] = {
      val arr = line.mkString.split("\t")
      val filter = arr(6)
      val ref = arr(3)
      val alt = arr(4)
      val coordination = arr(1).toInt
      val chromosome=arr(0)
      if (filter == "PASS" && ref.size == alt.size) Some((ref, coordination, alt,chromosome))
      else None
    }
  }

}
