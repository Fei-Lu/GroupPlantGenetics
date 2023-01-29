package daxing.v2.localAncestryInfer.simulation;

import java.time.LocalDate;
import java.util.List;

public class Order {
    private String orderNo;
    private LocalDate date;
    private String customerName;
    private List<OrderLine> orderLines;

    // Constructors, Getters, Setters and toString

    public Order(String orderNo, LocalDate localDate, String customerName, List<OrderLine> orderLines){
        this.orderNo = orderNo;
        this.date = localDate;
        this.customerName=customerName;
        this.orderLines=orderLines;
    }

    public Order(){

    }

    public String getOrderNo() {
        return orderNo;
    }

    public LocalDate getDate() {
        return date;
    }

    public String getCustomerName() {
        return customerName;
    }

    public List<OrderLine> getOrderLines() {
        return orderLines;
    }

    public void setCustomerName(String customerName) {
        this.customerName = customerName;
    }

    public void setDate(LocalDate date) {
        this.date = date;
    }

    public void setOrderNo(String orderNo) {
        this.orderNo = orderNo;
    }

    public void setOrderLines(List<OrderLine> orderLines) {
        this.orderLines = orderLines;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("Order{");
        sb.append("orderNo='").append(orderNo).append('\'');
        sb.append(", date=").append(date);
        sb.append(", customerName='").append(customerName).append('\'');
        sb.append(", orderLines=").append(orderLines);
        sb.append('}');
        return sb.toString();
    }
}
