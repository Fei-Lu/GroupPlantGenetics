package xiaohan.eQTL.javapractice;

public class basicdatatype {

    public basicdatatype(){
        this.ObjectWithdatatype();
        this.returnhashCode();
    }

    public void returnhashCode(){
        String sample = "whatthefuck";
        Object sample2 = sample.hashCode();
        System.out.println(sample2);
    }

    public void ObjectWithdatatype(){
        Integer a = 1;
        System.out.println(a+10);
    }

    public static void main(String[] args){
        new basicdatatype();
    }
}
