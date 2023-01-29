package daxing.v2.localAncestryInfer.simulation;

import java.math.BigDecimal;

public class OrderLine {
    private String item;
    private int quantity;
    private BigDecimal unitPrice;

    // Constructors, Getters, Setters and toString

    public OrderLine(String item, int quantity, BigDecimal unitPrice){
        this.item = item;
        this.quantity = quantity;
        this.unitPrice = unitPrice;
    }

    public OrderLine(){

    }

    public String getItem() {
        return item;
    }

    public int getQuantity() {
        return quantity;
    }

    public BigDecimal getUnitPrice() {
        return unitPrice;
    }

    public void setItem(String item) {
        this.item = item;
    }

    public void setQuantity(int quantity) {
        this.quantity = quantity;
    }

    public void setUnitPrice(BigDecimal unitPrice) {
        this.unitPrice = unitPrice;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("OrderLine{");
        sb.append("item='").append(item).append('\'');
        sb.append(", quantity=").append(quantity);
        sb.append(", unitPrice=").append(unitPrice);
        sb.append('}');
        return sb.toString();
    }
}
